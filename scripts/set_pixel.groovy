// Set pixel width to 0.25 microns for all images in project

def project = getProject()
def imageList = project.getImageList()

println "Setting pixel width to 0.25 microns for ${imageList.size()} images..."

for (entry in imageList) {
    def imageData = entry.readImageData()
    def server = imageData.getServer()
    
    // Set the image as current for the operations
    setBatchProjectAndImage(project, imageData)
    
    // Set pixel size to 0.25 microns
    setPixelSizeMicrons(0.25, 0.25)
    
    // Save
    entry.saveImageData(imageData)
    
    println "✓ ${entry.getImageName()}"
}

println "Done! All images now have 0.25 micron pixel width"