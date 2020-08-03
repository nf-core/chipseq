/*
 * Function to initialise default values and to generate a Groovy Map of module options
 */
def initOptions(Map args, String publish_dir) {
    def Map options = [:]
    options.args          = args.args ?: ''
    options.args2         = args.args2 ?: ''
    options.publish_by_id = args.publish_by_id ?: false
    options.publish_dir   = args.publish_dir ? args.publish_dir : publish_dir
    options.publish_files = args.publish_files ?: null
    options.suffix        = args.suffix ?: ''
    return options
}

/*
 * Function to save/publish module results
 *   if publish_files == null           : All files are published
 *   if publish_files == Map [:]        : No files are published
 *   if publish_files == Map [ext:path] : Only files that end with "ext" are published to "path" (appended to output directory)
 */
def saveFiles(filename, options, publish_dir) {
    if (!filename.endsWith('.version.txt')) {
        def publish_files = initOptions(options, publish_dir).publish_files
        if (publish_files instanceof Map) {
            for (ext in publish_files) {
                if (filename.endsWith(ext.key)) {
                    if (ext.value) {
                        return "${ext.value}/$filename"
                    } else {
                        return filename
                    }
                }
            }
        } else {
            return filename
        }
    }
}
