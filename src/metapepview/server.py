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

# In-memory session registry, match file to component_id
component_files = {}
upload_sessions = {}


def register_upload_endpoints(app):
    server = app.server

    @server.route("/upload-init", methods=["POST"])
    def upload_init():
        data = request.json
        filename = data["filename"]
        component_id = data["component_id"]

        # If new file replaces existing file from component, remove old file
        old_file = component_files.get(component_id)
        if old_file and os.path.exists(old_file):
            os.remove(old_file)

        upload_id = str(uuid.uuid4())
        temp_path = os.path.join(UPLOAD_FOLDER, upload_id + ".part")

        upload_sessions[upload_id] = {
            "filename": filename,
            "path": temp_path,
            "component_id": component_id

        }

        return jsonify({"upload_id": upload_id})


    @server.route("/upload-chunk", methods=["POST"])
    def upload_chunk():
        upload_id = request.form["upload_id"]
        chunk = request.files["chunk"]

        session = upload_sessions.get(upload_id)
        if not session:
            return "Invalid session", 400

        with open(session["path"], "ab") as f:
            f.write(chunk.read())

        return "OK", 200


    @server.route("/upload-complete", methods=["POST"])
    def upload_complete():
        data = request.json
        upload_id = data["upload_id"]

        session = upload_sessions.get(upload_id)
        if not session:
            return "Invalid session", 400
        
        component_id = session["component_id"]

        final_path = os.path.join(UPLOAD_FOLDER, session["filename"])
        os.rename(session["path"], final_path)

        component_files[component_id] = final_path

        del upload_sessions[upload_id]

        return jsonify({"filepath": final_path})


app = Dash(__name__, 
           external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
           suppress_callback_exceptions=True, 
           )

register_upload_endpoints(app)
