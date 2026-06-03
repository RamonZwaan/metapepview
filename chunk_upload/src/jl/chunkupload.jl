# AUTO GENERATED FILE - DO NOT EDIT

export chunkupload

"""
    chunkupload(;kwargs...)

A ChunkUpload component.

Keyword arguments:
- `id` (String; optional)
- `file_info` (Dict; optional)
- `progress` (Real; optional)
"""
function chunkupload(; kwargs...)
        available_props = Symbol[:id, :file_info, :progress]
        wild_props = Symbol[]
        return Component("chunkupload", "ChunkUpload", "chunk_upload", available_props, wild_props; kwargs...)
end

