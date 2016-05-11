# CVFinalProject
"Magic" frame project implementation for AIT-Budapest's Computer Vision Applications for Digital Cinema.  
TODO's marked in R2Image.cpp.
## Usage
src/imgpro input/0000000.jpg output/0000000.jpg -magic NUMBER_OF_IMAGES
## Step 1: Setup
### Iterate through images and call the correct functions on them
DONE: -magic tag in imgpro.cpp
NOTE: Coded based on the name format of jpg files provided by professor (e.g. 0000000.jpg, 0000001.jpg, etc).  
### Make tracking faster 
Something about using 4 threads? Might not be necessary...
## Step 2: Detect and store the "frame" boundary information
Probably need a separate function for this that will be used in steps 3 and 4  
## Step 3: Extract and store the frozen image from the first frame
magicExtractFrozen()
## Step 4: Replace the stuff inside the frame with frozen image 
DONE: magicReplaceFrameContent(nextImage)

## Methods and variables added for this project (incomplete):
### Methods:
**R2Image::magicFeature()** --> feature detection based on assignments from throughout the semester. Called on 
			a single image only.
			TODO add search for red/green/blue/black clusters around feature points to identify trackers
**R2Image::clusters(x,y)** --> given coordinates (x,y) of a feature point, check if surrounding pixels form something
			similar to the trackers we created. Return true if clusters of red/blue/green/black/white.
			TODO define thresholds for RGB values. 
**R2Image::findShiftedFrame(nextImage, prev_frame)** --> given frame coordinates prev_frame, returns the new coordinates
			of the frame in the nextImage, running a faster local search to locate the feature points
**R2Image::magicReplaceFrameContent(nextImage, shifted_frame)** --> takes in the next image and replaces the pixels 
			using the coordinates of the shifted frame calculated from findShiftedFrame
**R2Image::magicExtractFrozen()** --> Detects and stores the iformation of the first frame

### Members:
- Struct "Coordinates" to contain int values of x and y coordinates of pixels
- Array of 4 coordinates structs to contain pairs of ints (x,y), location of frame corners in each image.
Meant to be used as four center points for local searches when detecting the trackers in the next image.
- 2D array of doubles hMatrix, 3x3. Contains transformation between images. To be used for replacing image 
within the tracked frame. 

### Global vars.:
Empty R2Image freezeFrame. To be resized and set to the contents of the image within the bounds of the 
frame corners. Passed as parameter in replace image functions. Note: Image object freezeFrame must be a 
global variable, since Image class cannot contain itself.
(I don't think we need this according to Yi-tong's implementation? -Annie)