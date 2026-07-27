# AUTO GENERATED FILE - DO NOT EDIT

export chunkupload

"""
    chunkupload(;kwargs...)

A ChunkUpload component.

Keyword arguments:
- `id` (String; optional)
- `chunk_size` (Real; optional)
- `className` (String; optional)
- `file_info` (Dict; optional)
- `file_infos` (Array; optional)
- `multiple` (Bool; optional)
- `progress` (Real; optional)
- `style` (Dict; optional)
"""
function chunkupload(; kwargs...)
        available_props = Symbol[:id, :chunk_size, :className, :file_info, :file_infos, :multiple, :progress, :style]
        wild_props = Symbol[]
        return Component("chunkupload", "ChunkUpload", "chunk_upload", available_props, wild_props; kwargs...)
end

