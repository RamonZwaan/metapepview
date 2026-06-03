# AUTO GENERATED FILE - DO NOT EDIT

#' @export
chunkUpload <- function(id=NULL, file_info=NULL, progress=NULL) {
    
    props <- list(id=id, file_info=file_info, progress=progress)
    if (length(props) > 0) {
        props <- props[!vapply(props, is.null, logical(1))]
    }
    component <- list(
        props = props,
        type = 'ChunkUpload',
        namespace = 'chunk_upload',
        propNames = c('id', 'file_info', 'progress'),
        package = 'chunkUpload'
        )

    structure(component, class = c('dash_component', 'list'))
}
