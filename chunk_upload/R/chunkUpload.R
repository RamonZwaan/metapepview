# AUTO GENERATED FILE - DO NOT EDIT

#' @export
chunkUpload <- function(id=NULL, appearance=NULL, chunk_size=NULL, className=NULL, file_infos=NULL, multiple=NULL, progress=NULL, style=NULL, uploadText=NULL) {
    
    props <- list(id=id, appearance=appearance, chunk_size=chunk_size, className=className, file_infos=file_infos, multiple=multiple, progress=progress, style=style, uploadText=uploadText)
    if (length(props) > 0) {
        props <- props[!vapply(props, is.null, logical(1))]
    }
    component <- list(
        props = props,
        type = 'ChunkUpload',
        namespace = 'chunk_upload',
        propNames = c('id', 'appearance', 'chunk_size', 'className', 'file_infos', 'multiple', 'progress', 'style', 'uploadText'),
        package = 'chunkUpload'
        )

    structure(component, class = c('dash_component', 'list'))
}
