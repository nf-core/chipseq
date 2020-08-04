/*
 * Extract name of software from nf-core/modules process name using $task.process
 */
def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
}

/*
 * Function to initialise default values and to generate a Groovy Map of module options
 */
def initOptions(Map args) {
    def Map options = [:]
    options.args          = args.args ?: ''
    options.args2         = args.args2 ?: ''
    options.publish_by_id = args.publish_by_id ?: false
    options.publish_dir   = args.publish_dir ?: ''
    options.publish_files = args.publish_files ?: null
    options.suffix        = args.suffix ?: ''
    return options
}

/*
 * Tidy up and join elements of a list to return a path string
 */
def getPathFromList(path_list) {
  def paths = path_list.findAll { item -> !item?.trim().isEmpty() }  // Remove empty entries
  paths = paths.collect { it.trim().replaceAll("^[/]+|[/]+\$", "") } // Trim whitespace and trailing slashes
  return paths.join('/')
}

/*
 * Function to save/publish module results
 *   if publish_files == null           : All files are published
 *   if publish_files == Map [:]        : No files are published
 *   if publish_files == Map [ext:path] : Only files that end with "ext" are published to "path" (appended to output directory)
 */
def saveFiles(filename, options, publish_dir='', publish_id='') {
    if (!filename.endsWith('.version.txt')) {
        def ioptions = initOptions(options)
        def path_list = [ ioptions.publish_dir ?: publish_dir ]
        if (ioptions.publish_by_id) {
            path_list.add(publish_id)
        }
        if (ioptions.publish_files instanceof Map) {
            for (ext in ioptions.publish_files) {
                if (filename.endsWith(ext.key)) {
                    def ext_list = path_list.collect()
                    ext_list.add(ext.value)
                    return "/${getPathFromList(ext_list)}/$filename"
                }
            }
        } else {
            return "/${getPathFromList(path_list)}/$filename"
        }
    }
}
