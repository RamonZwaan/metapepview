import os
import uuid
from pathlib import Path
import shutil
from flask import request, jsonify
from dash import Dash
import dash_bootstrap_components as dbc
from metapepview.constants import GlobalConstants


# configure directory to store large uploads
UPLOAD_FOLDER = GlobalConstants.upload_folder
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

for item in Path(UPLOAD_FOLDER).iterdir():
    shutil.rmtree(item) if item.is_dir() else item.unlink()


def register_upload_endpoints(app):
    server = app.server

    # In-memory session registry, match component id to files
    server.component_files = {}

    # store file information during upload session
    server.upload_sessions = {}


    @server.route("/upload-init", methods=["POST"])
    def upload_init():
        """Trigger to start new upload session, delete old files related to
        component. 
        """
        data = request.json

        component_id = data["component_id"]

        # If new file replaces existing file from component, remove old file
        old_files = server.component_files.get(component_id)
        if old_files:
            for old_file in old_files:
                if os.path.exists(old_file):
                    os.remove(old_file)

        # set files connected to component id to empty
        server.component_files[component_id] = []

        return "OK", 200


    @server.route("/file-upload-init", methods=["POST"])
    def file_upload_init():
        """New file upload, create file entry in upload session.
        """
        data = request.json

        filename = data["filename"]
        component_id = data["component_id"]

        # create a new upload session on by-file basis
        upload_id = str(uuid.uuid4())
        temp_path = os.path.join(UPLOAD_FOLDER, upload_id + ".part")

        server.upload_sessions[upload_id] = {
            "filename": filename,
            "path": temp_path,
            "component_id": component_id
        }

        return jsonify({"upload_id": upload_id})


    @server.route("/upload-chunk", methods=["POST"])
    def upload_chunk():
        """Add chunk data to file
        """
        upload_id = request.form["upload_id"]
        chunk = request.files["chunk"]

        session = server.upload_sessions.get(upload_id)
        if not session:
            return "Invalid session", 400

        with open(session["path"], "ab") as f:
            f.write(chunk.read())

        return "OK", 200


    @server.route("/file-upload-complete", methods=["POST"])
    def file_upload_complete():
        """Upload for single file complete, store file in final path and close 
        
        """
        data = request.json
        upload_id = data["upload_id"]

        session = server.upload_sessions.get(upload_id)
        if not session:
            return "Invalid session", 400
        
        component_id = session["component_id"]

        final_path = os.path.join(UPLOAD_FOLDER, f"{component_id}_{session['filename']}")

        try:
            os.rename(session["path"], final_path)

            server.component_files[component_id].append(final_path)
        finally:
            del server.upload_sessions[upload_id]

        return jsonify({"path": final_path,
                        "filename": session['filename']})


app = Dash(__name__, 
           external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
           suppress_callback_exceptions=True, 
           )

# store large datasets server-side, storing key in dcc.Store
app.server.data_store = {}

register_upload_endpoints(app)
