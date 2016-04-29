# CVFinalProject
"Magic" frame project implementation for AIT-Budapest's Computer Vision Applications for Digital Cinema.
TODO's marked in R2Image.cpp.
## Usage
src/imgpro input/0000000.jpg output/0000000.jpg -magic NUMBER_OF_IMAGES
## Step 1: Setup
### Iterate through images and call the correct functions on them
DONE: READ COMMENTS IN imgpro.cpp (search for "magic" in the file).
NOTE: Coded based on the name format of jpg files provided by professor (e.g. 0000000.jpg, 0000001.jpg, etc).
### Make tracking faster (with 4 threads)
## Step 2: Detect and store the "frame" boundary information
Probably need a separate function for this that will be used in steps 3 and 4
## Step 3: Extract and store the frozen image from the first frame
magicExtractFrozen()
## Step 4: Replace the stuff inside the frame with frozen image 
magicReplaceFrameContent(nextImage)
