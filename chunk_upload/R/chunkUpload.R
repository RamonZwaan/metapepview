# AUTO GENERATED FILE - DO NOT EDIT

#' @export
chunkUpload <- function(id=NULL, chunk_size=NULL, className=NULL, file_info=NULL, file_infos=NULL, multiple=NULL, progress=NULL, style=NULL) {
    
    props <- list(id=id, chunk_size=chunk_size, className=className, file_info=file_info, file_infos=file_infos, multiple=multiple, progress=progress, style=style)
    if (length(props) > 0) {
        props <- props[!vapply(props, is.null, logical(1))]
    }
    component <- list(
        props = props,
        type = 'ChunkUpload',
        namespace = 'chunk_upload',
        propNames = c('id', 'chunk_size', 'className', 'file_info', 'file_infos', 'multiple', 'progress', 'style'),
        package = 'chunkUpload'
        )

    structure(component, class = c('dash_component', 'list'))
}
