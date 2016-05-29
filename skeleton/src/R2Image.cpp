// Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <iostream>
#include <time.h>
#include <math.h> //this is for ceil() and floor()



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////



R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
  // fit a 2D conic to five points
  R2Point p1(1.2,3.5);
  R2Point p2(2.1,2.2);
  R2Point p3(0.2,1.6);
  R2Point p4(0.0,0.5);
  R2Point p5(-0.2,4.2);

  // build the 5x6 matrix of equations
  double** linEquations = dmatrix(1,5,1,6);

  linEquations[1][1] = p1[0]*p1[0];
  linEquations[1][2] = p1[0]*p1[1];
  linEquations[1][3] = p1[1]*p1[1];
  linEquations[1][4] = p1[0];
  linEquations[1][5] = p1[1];
  linEquations[1][6] = 1.0;

  linEquations[2][1] = p2[0]*p2[0];
  linEquations[2][2] = p2[0]*p2[1];
  linEquations[2][3] = p2[1]*p2[1];
  linEquations[2][4] = p2[0];
  linEquations[2][5] = p2[1];
  linEquations[2][6] = 1.0;

  linEquations[3][1] = p3[0]*p3[0];
  linEquations[3][2] = p3[0]*p3[1];
  linEquations[3][3] = p3[1]*p3[1];
  linEquations[3][4] = p3[0];
  linEquations[3][5] = p3[1];
  linEquations[3][6] = 1.0;
  
  linEquations[4][1] = p4[0]*p4[0];
  linEquations[4][2] = p4[0]*p4[1];
  linEquations[4][3] = p4[1]*p4[1];
  linEquations[4][4] = p4[0];
  linEquations[4][5] = p4[1];
  linEquations[4][6] = 1.0;

  linEquations[5][1] = p5[0]*p5[0];
  linEquations[5][2] = p5[0]*p5[1];
  linEquations[5][3] = p5[1]*p5[1];
  linEquations[5][4] = p5[0];
  linEquations[5][5] = p5[1];
  linEquations[5][6] = 1.0;

  printf("\n Fitting a conic to five points:\n");
  printf("Point #1: %f,%f\n",p1[0],p1[1]);
  printf("Point #2: %f,%f\n",p2[0],p2[1]);
  printf("Point #3: %f,%f\n",p3[0],p3[1]);
  printf("Point #4: %f,%f\n",p4[0],p4[1]);
  printf("Point #5: %f,%f\n",p5[0],p5[1]);

  // compute the SVD
  double singularValues[7]; // 1..6
  double** nullspaceMatrix = dmatrix(1,6,1,6);
  svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

  // get the result
  printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],
    nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

  // make sure the solution is correct:
  printf("Equation #1 result: %f\n",  p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
                    p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
                    p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
                    p1[0]*nullspaceMatrix[4][smallestIndex] + 
                    p1[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #2 result: %f\n",  p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
                    p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
                    p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
                    p2[0]*nullspaceMatrix[4][smallestIndex] + 
                    p2[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #3 result: %f\n",  p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
                    p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
                    p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
                    p3[0]*nullspaceMatrix[4][smallestIndex] + 
                    p3[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #4 result: %f\n",  p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
                    p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
                    p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
                    p4[0]*nullspaceMatrix[4][smallestIndex] + 
                    p4[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #5 result: %f\n",  p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
                    p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
                    p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
                    p5[0]*nullspaceMatrix[4][smallestIndex] + 
                    p5[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  R2Point test_point(0.34,-2.8);

  printf("A point off the conic: %f\n", test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
                      test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
                      test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
                      test_point[0]*nullspaceMatrix[4][smallestIndex] + 
                      test_point[1]*nullspaceMatrix[5][smallestIndex] + 
                      nullspaceMatrix[6][smallestIndex]);

  return; 
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelHelper(double kernel[3][3])
{
  // Help from here: http://stackoverflow.com/questions/17815687/image-processing-implementing-sobel-filter
  R2Image temp(width, height);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      R2Pixel mag;
      for(int a = 0; a < 3; a++) {
        for(int b = 0; b < 3; b++) {
          int xn = i + a - 1;
          int yn = j + b - 1;
          mag += kernel[a][b] * Pixel(xn, yn);
        }
      }
      temp.Pixel(i,j) = mag;
    }
  }

  (*this) = temp;

}

void R2Image::
SobelX(void)
{
  // Apply the Sobel operator to the image in X direction

  double kernelX[3][3] = {{-1, 0, 1}, 
                          {-2, 0, 2}, 
                          {-1, 0, 1}};

  SobelHelper(kernelX);

}

void R2Image::
SobelY(void)
{
  // Apply the Sobel operator to the image in Y direction

  double kernelY[3][3] = {{-1, -2, -1}, 
                          {0, 0, 0}, 
                          {1, 2, 1}};

  SobelHelper(kernelY);
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}

// Non-linear filtering ////////////////////////////////////////////////
void R2Image::
Bilateral(double sigma)
{
  // Bilateral filter of the image
  printf("BILATERAL FILTER! OH!\n");
}

void R2Image::
Distort(void)
{
  // Apply Lens Distortion
  printf("LENS DISTORT! AH!\n");
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred
  
  R2Image temp(width, height);
  R2Image temp2(width, height);
  int kerSize = (int) (6 * sigma) + 1;
  double origin_dist = kerSize/2;
  double* sepKernel = new double[kerSize];

  // https://en.wikipedia.org/wiki/Gaussian_blur#Mechanics
  for (int i = 0; i < kerSize; i++){
    double fraction = 1 / (2 * 3.14 * pow(sigma, 2));
    double top_equ_part = pow((i - origin_dist), 2) + pow(origin_dist, 2);
    double exponent = -(top_equ_part) / (2 * pow(sigma, 2));
    sepKernel[i] = fraction * pow(2.718, exponent);
  }

  for (int y = origin_dist; y < height - origin_dist; y++) {
    for (int x = origin_dist; x < width - origin_dist; x++) {
      R2Pixel mag;
      double weightSum = 0.0;
      for (int j = 0; j < kerSize; j++){
        int yn = y + (j - origin_dist);
        mag += sepKernel[j] * Pixel(x, yn);
        weightSum += sepKernel[j];
      }
      temp.Pixel(x,y) = mag / weightSum;
    }
  }

  for (int y = origin_dist; y < height - origin_dist; y++) {
    for (int x = origin_dist; x < width - origin_dist; x++) {
      R2Pixel mag;
      double weightSum = 0.0;
      for (int i = 0; i < kerSize; i++){
        int xn = x + (i - origin_dist);
        mag += sepKernel[i] * temp.Pixel(xn, y);
        weightSum += sepKernel[i];
      }
      temp2.Pixel(x,y) = mag / weightSum;
    }
  }

  (*this) = temp2;
}

// CHANGED THIS TO R2IMAGE::POINT, CHANGED ALL OTHER METHODS ACCORDINGLY
// struct point {
//   int x;
//   int y;
//   double val;
// };


int compare(const void * p1, const void * p2) {
  if ((*(R2Image::point*)p1).val < (*(R2Image::point*)p2).val) {
    return 1;
  } else if ((*(R2Image::point*)p1).val > (*(R2Image::point*)p2).val) {
    return -1;
  } else {
    return 0;
  }
}

void R2Image::
Harris(double sigma)
{
  // TODO Optimize the features finder
  // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
  // Output should be 50% grey at flat regions, white at corners and black/dark near edges

  // Compute (I_x)^2, (I_y)^2, and I_x*I_y

  // timer goes here
  clock_t start, end;
  start = clock();

  R2Image Ix = R2Image(*this);
  Ix.SobelX();
  R2Image Iy = R2Image(*this);
  Iy.SobelY();
  R2Image Ix_Iy = R2Image(width, height);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Ix_Iy.Pixel(i,j) = Ix.Pixel(i,j) * Iy.Pixel(i,j);
      Ix.Pixel(i,j) = Ix.Pixel(i,j) * Ix.Pixel(i,j);
      Iy.Pixel(i,j) = Iy.Pixel(i,j) * Iy.Pixel(i,j);
    }
  }

  // Apply a Gaussian Blur of small kernel size to all three
  Ix.Blur(sigma);
  Iy.Blur(sigma);
  Ix_Iy.Blur(sigma);

  // Compute the Harris value for each pixel from these
  R2Image Harris = R2Image(width,height);
  // Store values in an array
  point *values = new point[width * height];
  // For normalization such that 0 becomes grey
  R2Pixel grey = R2Pixel(0.5, 0.5, 0.5, 0.0);

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      R2Pixel det = Ix.Pixel(i,j) * Iy.Pixel(i,j) - Ix_Iy.Pixel(i,j) * Ix_Iy.Pixel(i,j);
      R2Pixel trace = Ix.Pixel(i,j) + Iy.Pixel(i,j);
      double alpha = 0.04;
      Harris.Pixel(i,j) = det - alpha * (trace * trace);
      Harris.Pixel(i,j) += grey;
      Harris.Pixel(i,j);
      point curPoint;
      curPoint.val = Harris.Pixel(i,j).Red() + Harris.Pixel(i,j).Green() + Harris.Pixel(i,j).Blue();
      curPoint.x = i;
      curPoint.y = j;
      values[i * height + j] = curPoint;
    }
  }

  // Sort the array according to values
  qsort(values, width * height, sizeof(point), compare);

  // Find 150 features with very high corner score
  point *topValues = new point[150];
  int featureCount = 0;

  for (int idx = 0; featureCount < 150 && idx < width*height; idx++) {
    point curPoint = values[idx];
    bool close = false;
    // Make sure there is at least 10 pixel distance between them
    for (int feat = 0; feat < featureCount; feat++) {
      int x = topValues[feat].x - curPoint.x;
      int y = topValues[feat].y - curPoint.y;
      if (sqrt(x*x + y*y) < width/50) {
        close = true;
      }
    }
    if(close) {
      continue;
    }
    topValues[featureCount++] = curPoint;
  }

  end = clock();
  float diff = (end - start) / CLOCKS_PER_SEC;
  std::cout << "Run time for Annie's is " << diff << std::endl;

  // Mark these points on the output image (red cross)
  R2Pixel red = R2Pixel(1.0, 0.0, 0.0, 0.0);
  for (int i = 0; i < 150; i++) {
    point curPoint = topValues[i];
    Pixel(curPoint.x, curPoint.y) = red;
    for (int i = 0; i < 5; i++) {
      Pixel(curPoint.x - i, curPoint.y) = red;
      Pixel(curPoint.x, curPoint.y  - i) = red;
      Pixel(curPoint.x + i, curPoint.y) = red;
      Pixel(curPoint.x, curPoint.y + i) = red;
    }
  }

  (*this) = Harris;
}


void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.

  // Gaussian is a low-pass filter
  R2Image original(width, height);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      original.Pixel(i,j) = Pixel(i,j);
    }
  }

  // Make this blurred
  Blur(2);

  R2Image high(width, height);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      high.Pixel(i,j) = original.Pixel(i,j) - Pixel(i,j);
    }
  }

  R2Image final(width, height);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      final.Pixel(i,j) = high.Pixel(i,j) * 2 + Pixel(i,j);
      final.Pixel(i,j).Clamp();
    }
  }

  (*this) = final;

}

void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if(x0>x1) {
    int x=y1;
    y1=y0;
    y0=x;
    x=x1;
    x1=x0;
    x0=x;
  }
  int deltax = x1 - x0;
  int deltay = y1 - y0;
  float error = 0;
  float deltaerr = 0.0;
  if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));
  int y = y0;
  for(int x=x0;x<=x1;x++) {
    Pixel(x,y).Reset(r,g,b,1.0);
    error = error + deltaerr;
    if(error>=0.5) {
      if(deltay>0) y = y + 1;
      else y = y - 1;
      error = error - 1.0;
    }
  }
  if(x0>3 && x0<width-3 && y0>3 && y0<height-3) {
    for(int x=x0-3;x<=x0+3;x++){
        for(int y=y0-3;y<=y0+3;y++){
          Pixel(x,y).Reset(r,g,b,1.0);
        }
    }
  }
}

struct translation {
  int xOriginal;
  int yOriginal;
  int xTranslated;
  int yTranslated;
};


double** hMatrixCalculation(translation translation[4]) {
  // Build the 8*9 matrix of equations for matrix A
  double** matrixP = dmatrix(1,8,1,9);
  // Loop over all point pairs and add two rows to matrixP
  for (int i = 1; i < 5; i++) {
      int x_a = translation[i-1].xOriginal;
      int y_a = translation[i-1].yOriginal;
      int w_a = 1;
      int x_b = translation[i-1].xTranslated;  //x'
      int y_b = translation[i-1].yTranslated;  //y'
      int w_b = 1;                             //w'

      // [0 0 0 -w_b*x_a -w_b*y_a -w_b*w_a y_b*x_a y_b*y_a y_b*w_a]
      matrixP[i*2-1][1] = 0; // 1
      matrixP[i*2-1][2] = 0; // 2
      matrixP[i*2-1][3] = 0; // 3
      matrixP[i*2-1][4] = -w_b*x_a; // 4
      matrixP[i*2-1][5] = -w_b*y_a; // 5
      matrixP[i*2-1][6] = -w_b*w_a; // 6
      matrixP[i*2-1][7] = y_b*x_a; // 7
      matrixP[i*2-1][8] = y_b*y_a; // 8
      matrixP[i*2-1][9] = y_b*w_a; // 9
      // [w_b*x_a w_b*y_a w_b*w_a 0 0 0 -x_b*x_a -x_b*y_a -x_b*w_a]
      matrixP[i*2][1] = w_b*x_a; // 1
      matrixP[i*2][2] = w_b*y_a; // 2
      matrixP[i*2][3] = 1;//w_b*w_a; // 3
      matrixP[i*2][4] = 0; // 4
      matrixP[i*2][5] = 0; // 5
      matrixP[i*2][6] = 0; // 6
      matrixP[i*2][7] = -x_b*x_a; // 7
      matrixP[i*2][8] = -x_b*y_a; // 8
      matrixP[i*2][9] = -x_b*w_a; // 9
  }
  // Find the nullspace of the matrix A (see the SVD example in the code)
  // Compute the SVD
  double singularValues[10];
  double** nullSpaceMatrix = dmatrix(1,9,1,9);
  svdcmp(matrixP, 8, 9, singularValues, nullSpaceMatrix);
  // Find the smallest singular value:
  int smallestIndex = 1;

  for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;
  // printf("SINGULAR VALUES: ");
  //   for(int i=1;i<10;i++) printf(" %f ", singularValues[i]);
  // printf("\n");
    
  // Build the matrix H_norm from the resulting 9 numbers
  double** H = dmatrix(1,3, 1, 3);
  H[1][1] = nullSpaceMatrix[1][smallestIndex];
  H[1][2] = nullSpaceMatrix[2][smallestIndex];
  H[1][3] = nullSpaceMatrix[3][smallestIndex];
  H[2][1] = nullSpaceMatrix[4][smallestIndex];
  H[2][2] = nullSpaceMatrix[5][smallestIndex];
  H[2][3] = nullSpaceMatrix[6][smallestIndex];
  H[3][1] = nullSpaceMatrix[7][smallestIndex];
  H[3][2] = nullSpaceMatrix[8][smallestIndex];
  H[3][3] = nullSpaceMatrix[9][smallestIndex];

  return H;
}

void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
  R2Image orig = R2Image(*this);
  // Detect 150 features on image A (done)
  // Compute (I_x)^2, (I_y)^2, and I_x*I_y
  int sigma = 3;
  R2Image Ix = R2Image(*this);
  Ix.SobelX();
  R2Image Iy = R2Image(*this);
  Iy.SobelY();
  R2Image Ix_Iy = R2Image(width, height);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Ix_Iy.Pixel(i,j) = Ix.Pixel(i,j) * Iy.Pixel(i,j);
      Ix.Pixel(i,j) = Ix.Pixel(i,j) * Ix.Pixel(i,j);
      Iy.Pixel(i,j) = Iy.Pixel(i,j) * Iy.Pixel(i,j);
    }
  }
  // Apply a Gaussian Blur of small kernel size to all three
  Ix.Blur(sigma);
  Iy.Blur(sigma);
  Ix_Iy.Blur(sigma);
  // Compute the Harris value for each pixel from these
  R2Image Harris = R2Image(width,height);
  // Store values in an array
  point *values = new point[width * height];
  // For normalization such that 0 becomes grey
  R2Pixel grey = R2Pixel(0.5, 0.5, 0.5, 0.0);
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      R2Pixel det = Ix.Pixel(i,j) * Iy.Pixel(i,j) - Ix_Iy.Pixel(i,j) * Ix_Iy.Pixel(i,j);
      R2Pixel trace = Ix.Pixel(i,j) + Iy.Pixel(i,j);
      double alpha = 0.04;
      Harris.Pixel(i,j) = det - alpha * (trace * trace);
      Harris.Pixel(i,j) += grey;
      Harris.Pixel(i,j);
      point curPoint;
      curPoint.val = Harris.Pixel(i,j).Red() + Harris.Pixel(i,j).Green() + Harris.Pixel(i,j).Blue();
      curPoint.x = i;
      curPoint.y = j;
      values[i * height + j] = curPoint;
    }
  }
  // Sort the array according to values
  qsort(values, width * height, sizeof(point), compare);
  // Find 150 features with very high corner score
  printf("FINDING FEATURES\n");
  point *topValues = new point[150];
  int featureCount = 0;
  for (int idx = 0; featureCount < 150 && idx < width*height; idx++) {
    point curPoint = values[idx];
    bool close = false;
    // Make sure there is at least 10 pixel distance between them
    for (int feat = 0; feat < featureCount; feat++) {
      int x = topValues[feat].x - curPoint.x;
      int y = topValues[feat].y - curPoint.y;
      if (sqrt(x*x + y*y) < 10) {
        close = true;
      }
    }
    if(close) {
      continue;
    }
    topValues[featureCount++] = curPoint;
  }

   R2Pixel BLUE = R2Pixel(0,0,1,1);
  for (int i = 139; i < 150; i++) {
    point cur = topValues[i];
    for (int j = -5; j < 6; j++) {
      for(int z = -5; z < 6; z++){
        Pixel(cur.x + j, cur.y + z) = BLUE;
      }
    }
    
  }

  // // Search area size should be reasonable (20% image size)
  // printf("FINDING TRANSLATIONS\n");
  // translation *translations = new translation[150];
  // int searchWidth = width/5;
  // int searchHeight = height/5;
  // const int offset = 13;
  // // For each point run a local search
  // for (int i = 0; i < 150; i++) {
  //   point curPoint = topValues[i];
  //   // Search loop should check SSD at each location within the region
  //   int minX = curPoint.x;
  //   int minY = curPoint.y;
  //   double minSum = 999999.0;
  //   for (int x = curPoint.x - searchWidth; x < curPoint.x + searchWidth; x++) {
  //     for (int y = curPoint.y - searchHeight; y < curPoint.y + searchHeight; y++) {
  //       double sum = 0.0;
  //       for (int i = -offset/2; i < offset/2; i++) {
  //         for (int j = -offset/2; j < offset/2; j++) {
  //           // check for bounds
  //           int xpp = curPoint.x + i;
  //           int ypp = curPoint.y + j;
  //           if (xpp < 0 || xpp > width) {
  //             xpp = curPoint.x;
  //           }
  //           if (ypp < 0 || ypp > height) {
  //             ypp = curPoint.y;
  //           }
  //           int xp = x + i;
  //           int yp = y + j;
  //           if (xp < 0 || xp > width) {
  //             xp = curPoint.x;
  //           }
  //           if (yp < 0 || yp > height) {
  //             yp = curPoint.y;
  //           }
  //           R2Pixel diff = orig.Pixel(xpp, ypp) - otherImage->Pixel(xp, yp);
  //           sum += diff.Red()*diff.Red() + diff.Green()*diff.Green() + diff.Blue()*diff.Blue();
  //         }
  //       }
  //       if (sum < minSum) {
  //         minX = x;
  //         minY = y;
  //         minSum = sum;
  //       }
  //     }
  //   }
  //   // Store the translation vector
  //   translation curTrans;
  //   curTrans.xOriginal = curPoint.x;
  //   curTrans.yOriginal = curPoint.y;
  //   curTrans.xTranslated = minX;
  //   curTrans.yTranslated = minY;
  //   translations[i] = curTrans;
  // }

  // // RANSAC
  // printf("RANSAC\n");
  // const double threshold = 4;
  // int ransacTrials = 200;
  // double** bestH = dmatrix(1,3, 1, 3);
  // int mostInliers = 0;
  // while(ransacTrials > 0) {
  //   // Randomly select a small subset of points (at least 4)
  //   translation *fourRandPoints = new translation[4];
  //   int randNums[150];
  //   for (int i = 0; i < 150; i++) {
  //     randNums[i] = i;
  //   }
  //   std::random_shuffle(std::begin(randNums), std::end(randNums));
  //   for (int i = 0; i < 4; i++) {
  //     fourRandPoints[i] = translations[randNums[i]];
  //   }
  //   double** hNorm = hMatrixCalculation(fourRandPoints);
  //   double** H = hNorm;
  //   int inliers = 0;
  //   for (int i = 0; i < 150; i++) {
  //     double final[3];
  //     for (int j = 1; j < 4; ++j) {
  //         final[j-1] = H[j][1] * translations[i].xOriginal + H[j][2] * translations[i].yOriginal + H[j][3] * 1;
  //     }
  //     // Project points x to x' to find matches
  //     double hypothesisX = final[0]/final[2];
  //     double hypothesisY = final[1]/final[2];
  //     double xComponent = translations[i].xTranslated - hypothesisX;
  //     double yComponent = translations[i].yTranslated - hypothesisY;
  //     double diff = sqrt(xComponent*xComponent + yComponent*yComponent);
  //     if (diff < threshold) {
  //       inliers++;
  //     }
  //   }
  //   if(inliers > mostInliers) {
  //     bestH = H;
  //     mostInliers = inliers;
  //   }
  //   ransacTrials -= 1;
  // }

  // // After N trials, choose the largest inlier set and re-estimate the translation
  // printf("COLOR TRANSLATIONS\n");
  // for (int i = 0; i < 150; i++) {
  //   // Determine inliers
  //   translation curTrans = translations[i];
  //   double final[3];
  //   for (int i = 1; i < 4; ++i) {
  //       final[i-1] = bestH[i][1] * curTrans.xOriginal + bestH[i][2] * curTrans.yOriginal + bestH[i][3] * 1;
  //   }
  //   double hypothesisX = final[0]/final[2];
  //   double hypothesisY = final[1]/final[2];
  //   double xComponent = curTrans.xTranslated - hypothesisX;
  //   double yComponent = curTrans.yTranslated - hypothesisY;
  //   double diff = sqrt(xComponent*xComponent + yComponent*yComponent);
  //   if (diff < threshold) {
  //     // inlier is green
  //     line((int)curTrans.xOriginal, (int)curTrans.xTranslated, (int)curTrans.yOriginal, (int)curTrans.yTranslated, 0.0, 1.0, 0.0);
  //   } else {
  //     // outlier is red
  //     line((int)curTrans.xOriginal, (int)curTrans.xTranslated, (int)curTrans.yOriginal, (int)curTrans.yTranslated, 1.0, 0.0, 0.0);
  //   }
  // }
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
  // Homography estimation phase 1
  // Implement a function that computes the H Homography matrix from these point correspondences:
  // A:
  //     (0,0)->(1,0)
  //     (1,0)->(2,0)
  //     (1,1)->(2,1)
  //     (0,1)->(1,1)
  // B:
  //     (0,0)->(1,2)
  //     (1,0)->(1,1)
  //     (1,1)->(3,1)
  //     (0,1)->(3,2)
  translation *transA = new translation[4];
  translation A1;
  A1.xOriginal = 0;
  A1.yOriginal = 0;
  A1.xTranslated = 1;
  A1.yTranslated = 0;
  transA[0] = A1;
  translation A2;
  A2.xOriginal = 1;
  A2.yOriginal = 0;
  A2.xTranslated = 2;
  A2.yTranslated = 0;
  transA[1] = A2;
  translation A3;
  A3.xOriginal = 1;
  A3.yOriginal = 1;
  A3.xTranslated = 2;
  A3.yTranslated = 1;
  transA[2] = A3;
  translation A4;
  A4.xOriginal = 0;
  A4.yOriginal = 1;
  A4.xTranslated = 1;
  A4.yTranslated = 1;
  transA[3] = A4;
  
  translation *transB = new translation[4];
  translation B1;
  B1.xOriginal = 0;
  B1.yOriginal = 0;
  B1.xTranslated = 1;
  B1.yTranslated = 2;
  transB[0] = B1;
  translation B2;
  B2.xOriginal = 1;
  B2.yOriginal = 0;
  B2.xTranslated = 1;
  B2.yTranslated = 1;
  transB[1] = B2;
  translation B3;
  B3.xOriginal = 1;
  B3.yOriginal = 1;
  B3.xTranslated = 3;
  B3.yTranslated = 1;
  transB[2] = B3;
  translation B4;
  B4.xOriginal = 0;
  B4.yOriginal = 1;
  B4.xTranslated = 3;
  B4.yTranslated = 2;
  transB[3] = B4;

  double** HA = hMatrixCalculation(transA);
  double** HB = hMatrixCalculation(transB);

  // Print solution
  printf("H MATRIX: A\n");
  for (int i = 1; i < 4; ++i) {
    for (int j = 1; j < 4; ++j) {
      std::cout << HA[i][j] << ' ';
    }
    std::cout << std::endl;
  }

  // Check solution
  for (int i = 0; i < 4; ++i) {
    translation trans = transA[i];
    int x = trans.xOriginal;
    int y = trans.yOriginal;
    printf("ORIGINAL (%d, %d) * H should be (%d, %d)\n", x, y, trans.xTranslated, trans.yTranslated);
    double final[3];
    for (int i = 1; i < 4; ++i) {
        final[i-1] = HA[i][1] * x + HA[i][2] * y + HA[i][3] * 1;
    }
    printf("Correct? (%f, %f)\n", final[0]/final[2], final[1]/final[2]);
  }

  printf("H MATRIX: B\n");
  for (int i = 1; i < 4; ++i) {
    for (int j = 1; j < 4; ++j) {
      std::cout << HB[i][j] << ' ';
    }
    std::cout << std::endl;
  }

  for (int i = 0; i < 4; ++i) {
    translation trans = transB[i];
    int x = trans.xOriginal;
    int y = trans.yOriginal;
    printf("ORIGINAL (%d, %d) * H should be (%d, %d)\n", x, y, trans.xTranslated, trans.yTranslated);
    double final[3];
    for (int i = 1; i < 4; ++i) {
        final[i-1] = HB[i][1] * x + HB[i][2] * y + HB[i][3] * 1;
    }
    printf("Correct? (%f, %f)\n", final[0]/final[2], final[1]/final[2]);
  }
}





/////////////////////////
/* FRAME FINAL PROJECT */
/////////////////////////




R2Image freezeFrame; 

// constants for identifying trackers
double threshold = 0.3; // color -- need to tighten this
double blueThreshold = 0.2;
double greenThreshold = 0.04;
int radius = 3; // pixels to search around the center -- potentially widen this

double sigma = 3.0;
double alpha = 0.04;
int numFeaturePoints = 300;
int minFeatureDistance = 10;

R2Pixel red = R2Pixel(1,0,0,1);
R2Pixel green = R2Pixel(0,1,0,1);
R2Pixel blue = R2Pixel(0,0,1,1);;

vector<R2Image::coordinates> redTracker;
vector<R2Image::coordinates> greenTracker;
vector<R2Image::coordinates> blueTracker;

struct compare {
  bool operator() (const R2Image::coordinates &c1, const R2Image::coordinates &c2) { return (c1.y < c2.y); }
} comparison;



void R2Image::magicFeature(void) {

  /* HARRIS FEATURE DETECTOR */ 

  // Compute (I_x)^2, (I_y)^2, and I_x*I_y
  R2Image Ix = R2Image(*this);
  Ix.SobelX();
  R2Image Iy = R2Image(*this);
  Iy.SobelY();

  R2Image Ix_Iy = R2Image(width, height);

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Ix_Iy.Pixel(i,j) = Ix.Pixel(i,j) * Iy.Pixel(i,j);
      Ix.Pixel(i,j) = Ix.Pixel(i,j) * Ix.Pixel(i,j);
      Iy.Pixel(i,j) = Iy.Pixel(i,j) * Iy.Pixel(i,j);
    }
  }

  // Apply a Gaussian Blur of small kernel size to all three
  Ix.Blur(sigma);
  Iy.Blur(sigma);
  Ix_Iy.Blur(sigma);

  // Compute the Harris value for each pixel from these
  R2Image Harris = R2Image(width,height);
  
  // Store values in an array
  point *values = new point[width * height];
  
  // For normalization such that 0 becomes grey
  R2Pixel grey = R2Pixel(0.5, 0.5, 0.5, 0.0);

  for (int i = 0; i < width; i++) {

    for (int j = 0;  j < height; j++) {
    
      R2Pixel det = Ix.Pixel(i,j) * Iy.Pixel(i,j) - Ix_Iy.Pixel(i,j) * Ix_Iy.Pixel(i,j);
      R2Pixel trace = Ix.Pixel(i,j) + Iy.Pixel(i,j);
      Harris.Pixel(i,j) = det - alpha * (trace * trace);
      Harris.Pixel(i,j) += grey;
      Harris.Pixel(i,j);
      point curPoint;
      curPoint.val = Harris.Pixel(i,j).Red() + Harris.Pixel(i,j).Green() + Harris.Pixel(i,j).Blue();
      curPoint.x = i;
      curPoint.y = j;
      values[i * height + j] = curPoint;
    
    }
  }

  // Sort the array according to values
  qsort(values, width * height, sizeof(point), compare);
  
  int featureCount = 0;
  for (int idx = 0; featureCount < numFeaturePoints && idx < width*height; idx++) {
    point curPoint = values[idx];
    bool close = false;
    // Make sure there is at least 10 pixel distance between them
    for (int feat = 0; feat < featureCount; feat++) {
  
      int x = this->topValues[feat].x - curPoint.x;
      int y = this->topValues[feat].y - curPoint.y;
  
      if (sqrt(x*x + y*y) < minFeatureDistance) {
        close = true;
      }
  
    }
  
    if(close) {
      continue;
  
    }
  
    this->topValues[featureCount++] = curPoint;
  }


  R2Image::coordinates temp = {};

  int num_foundPoints = 0;
  for (int i = 0; i < numFeaturePoints; i++) {
    temp.x = this->topValues[i].x;
    temp.y = this->topValues[i].y;
    
    if (clusters(temp, red)) {
      num_foundPoints += 1;
      redTracker.push_back(temp);
      continue;
    }

    if (clusters(temp, green)) {
      num_foundPoints += 1;
      greenTracker.push_back(temp);
      continue;
    }

    if (clusters(temp, blue)) {
      num_foundPoints += 1;
      blueTracker.push_back(temp);

      continue;
    }

  }

  sort(redTracker.begin(), redTracker.end(), comparison);

  int last_big_jump = 1; //stores the biggest jump. between redTracker[last_big_jump] and redTracker[last_big_jump - 1]
  for (int i = 2; i < redTracker.size(); i++) {
    if((redTracker[last_big_jump].y - redTracker[last_big_jump - 1].y) < (redTracker[i].y - redTracker[i - 1].y)){
      last_big_jump = i;
    }
  }

  double blueX = 0;
  double blueY = 0;
  for (int i = 0; i < blueTracker.size(); i++) {
      blueX += blueTracker[i].x;
      blueY += blueTracker[i].y;
    }

  blueX /= blueTracker.size();
  blueY /= blueTracker.size();

  this->frame_corners[0].x = static_cast<int> (blueX);
  this->frame_corners[0].y = static_cast<int> (blueY);

  /* sort values stored in red tracker by y-value;
  look for the biggest jump; the points before the biggest jump are
  in the lower red tracker (bottom right corner), frame corners index
  value is 2 */

  double redX1 = 0;
  double redY1 = 0;
  for (int i = last_big_jump; i < redTracker.size(); i++) {
    redX1 += redTracker[i].x;
    redY1 += redTracker[i].y;
  }

  redX1 /= (redTracker.size() - last_big_jump);
  redY1 /= (redTracker.size() - last_big_jump);

  this->frame_corners[1].x = static_cast<int> (redX1);
  this->frame_corners[1].y = static_cast<int> (redY1);

  double redX2 = 0;
  double redY2 = 0;
  for (int i = 0; i < last_big_jump; i++) {
    redX2 += redTracker[i].x;
    redY2 += redTracker[i].y;
  }

  redX2 /= last_big_jump;
  redY2 /= last_big_jump;

  this->frame_corners[2].x = static_cast<int> (redX2);
  this->frame_corners[2].y = static_cast<int> (redY2);

  double greenX = 0;
  double greenY = 0;
  for (int i = 0; i < greenTracker.size(); i++) {
      greenX += greenTracker[i].x;
      greenY += greenTracker[i].y;
    }

  greenX /= greenTracker.size();
  greenY /= greenTracker.size();

  this->frame_corners[3].x = static_cast<int> (greenX);
  this->frame_corners[3].y = static_cast<int> (greenY);

  redTracker.clear();
  greenTracker.clear();
  blueTracker.clear();

}

bool R2Image::clusters(coordinates center, R2Pixel color) {

  int num_redPoints = 0;
  int num_greenPoints = 0;
  int num_bluePoints = 0;

  // detect search color
  bool red = false;
  bool green = false;
  bool blue = false;

  if (color.Red() > 0) {
    red = true;
  } else if (color.Green() > 0) {
    green = true;
  } else if (color.Blue() > 0) {
    blue = true;
  } 

  for (int i = -1 * radius; i < radius + 1; i++) {
    for (int j = -1 * radius; j < radius + 1; j++) {

      // make sure search doesn't fall of the edge of the image
      if (i + center.x < 0 || i + center.x > width || j + center.y < 0 || j + center.y > height) {
        continue;
      }

      // look for red clusters
      if (red &&
        Pixel(i + center.x,j + center.y).Red() - Pixel(i + center.x,j + center.y).Blue() >= threshold &&
        Pixel(i + center.x,j + center.y).Red() - Pixel(i + center.x,j + center.y).Green() >= threshold) {
        
        num_redPoints += 1;
        continue;
      }

      // look for green clusters
      if (green &&
        Pixel(i + center.x,j + center.y).Green() - Pixel(i + center.x,j + center.y).Blue() >= greenThreshold &&
        Pixel(i + center.x,j + center.y).Green() - Pixel(i + center.x,j + center.y).Red() >= greenThreshold) {
        num_greenPoints += 1;
        continue;
      }

      // look for blue clusters
      if (blue &&
        Pixel(i + center.x,j + center.y).Blue() - Pixel(i + center.x,j + center.y).Red() >= blueThreshold &&
        Pixel(i + center.x,j + center.y).Blue() - Pixel(i + center.x,j + center.y).Green() >= blueThreshold) {
        
        num_bluePoints += 1;
        continue;
      }

    }

  }


  if (num_redPoints > 0 || num_greenPoints > 0 || num_bluePoints > 0) {
    return true;
  }

  return false;
}

R2Image::frame R2Image::
findShiftedFrame(R2Image * prevImage, R2Image * nextImage, coordinates prev_frame[4]){
  // Will run a local search using the previous frame coordinated to find new shifted frame
  frame new_frame;

  // Search area size should be reasonable
  int searchWidth = width/5;
  int searchHeight = height/5;
  const int offset = 17;
  // For each point run a local search with SSD
  for(int a = 0; a < 4; a++){
    coordinates center;
    coordinates prevCoord = prev_frame[a];
    int minX = prevCoord.x;
    int minY = prevCoord.y;
    double minSum = 999999.0;

    int prev_x = prevCoord.x;
    int prev_y = prevCoord.y;
    int cur_x = prevCoord.x;
    int cur_y = prevCoord.y;

    for (int x = prevCoord.x - searchWidth; x < prevCoord.x + searchWidth; x++) {
      for (int y = prevCoord.y - searchHeight; y < prevCoord.y + searchHeight; y++) {
        double sum = 0.0;
        for (int i = -offset/2; i < offset/2; i++) {
          for (int j = -offset/2; j < offset/2; j++) {
            // Check for bounds
            if (prevCoord.x + i > 0 && prevCoord.x + i < width) {prev_x = prevCoord.x + i;}
            if (prevCoord.y + j > 0 && prevCoord.y + j < height) {prev_y = prevCoord.y + j;}
            if (x + i > 0 && x + i < width) {cur_x = x + i;}
            if (y + j > 0 && y + j < height) {cur_y = y + j;}
            R2Pixel diff = prevImage->Pixel(prev_x, prev_y) - nextImage->Pixel(cur_x, cur_y);
            sum += diff.Red()*diff.Red() + diff.Green()*diff.Green() + diff.Blue()*diff.Blue();
          }
        }
        if (sum < minSum) {
          minX = x;
          minY = y;
          minSum = sum;
        }
      }
    }
    // Store new point here
    center.x = minX;
    center.y = minY;
    new_frame.coordinates[a] = center;
  }

  return new_frame;
}



R2Image::line_equation computeLine(R2Image::coordinates a, R2Image::coordinates b){
  R2Image::line_equation answer = {}; //initalizing all values in answer to zero.
  //important so that answer.vertical is initialized to 0.
  if(a.x == b.x){
    answer.vertical = true;
    answer.b = a.x;
  }else{
    answer.m = (a.y - b.y)/(double)(a.x - b.x);
    answer.b = a.y - a.x * answer.m;
  }
  return answer;
}

R2Image::d_coordinates matrixMult(double x_val, double y_val, double** H){
  R2Image::d_coordinates answer = {};
  double x_prime_val;
  double y_prime_val;
  double z_prime_val;
 
  x_prime_val  = H[1][1] * x_val;
  x_prime_val += H[1][2] * y_val;
  x_prime_val += H[1][3];

  y_prime_val  = H[2][1] * x_val;
  y_prime_val += H[2][2] * y_val;
  y_prime_val += H[2][3];

  z_prime_val  = H[3][1] * x_val;
  z_prime_val += H[3][2] * y_val;
  z_prime_val += H[3][3];

//  fprintf(stdout, "z_prime_val %f for x: %d and y: %d  \n", z_prime_val, x_val, y_val);
  answer.x = x_prime_val / z_prime_val;
  answer.y = y_prime_val / z_prime_val;
  return answer;
}

void R2Image::
magicReplaceFrameContent(R2Image * nextImage, frame sh_frame)
{
  translation * trans = new translation[4]; //this'll store the upper left, upper right,lower right, lower left
  //translations in that order (also if you think of a more descriptive name than "trans" I'm all ears.)
  

  //shifted_frame should store the points of the frame after shifting.

  coordinates shifted_frame[4];

  shifted_frame[0].x = sh_frame.coordinates[0].x;
  shifted_frame[0].y = sh_frame.coordinates[0].y;

  shifted_frame[1].x = sh_frame.coordinates[1].x;
  shifted_frame[1].y = sh_frame.coordinates[1].y;

  shifted_frame[2].x = sh_frame.coordinates[2].x;
  shifted_frame[2].y = sh_frame.coordinates[2].y;
  
  shifted_frame[3].x = sh_frame.coordinates[3].x;
  shifted_frame[3].y = sh_frame.coordinates[3].y;

 /*
  0-------1
  |       | 
  |       |
  |       |
  3-------2
  */

  //sneaky sneaky. In order to get the matrix that warps nextImage to original image I plug the numbers in backward.
  for(int i = 0; i < 4; i++){
    //frame_corners stores the upper left, upper right, lower left, lower right positions of the frame corners.
    trans[i].xTranslated = frame_corners[i].x;
    trans[i].yTranslated = frame_corners[i].y;
    trans[i].xOriginal = shifted_frame[i].x;
    trans[i].yOriginal = shifted_frame[i].y;
  }

  // for(int i = 0; i< 4; i++){
  //   fprintf(stdout, "trans[%d] original  x: %d y: %d \n      
        //translated x: %d y %d \n", i,trans[i].xOriginal, trans[i].yOriginal, trans[i].xTranslated, trans[i].yTranslated);
  // }

  double** H = hMatrixCalculation(trans);
  //important to note: This H should map points from the nextImage to the original image
  //not the other way round. This is so we can inverse warp.


  //all bets are off for the orientation of the corners.
  //the frame could have been twirled around a couple times.
  int leftmost = shifted_frame[0].x;
  int rightmost = shifted_frame[0].x;
  int topmost = shifted_frame[0].y;
  int botmost = shifted_frame[0].y;
  for(int i = 1; i < 4; i++){
    if(shifted_frame[i].x < leftmost){
      leftmost = shifted_frame[i].x;
    }
    if(shifted_frame[i].x > rightmost){
      rightmost = shifted_frame[i].x;
    }
    if(shifted_frame[i].y > topmost){
      topmost = shifted_frame[i].y;
    }
    if(shifted_frame[i].y < botmost){
      botmost = shifted_frame[i].y;
    }
  }

  R2Image::line_equation frame[4];
  for(int i = 0; i < 4; i++){
    frame[i] = computeLine(shifted_frame[i], shifted_frame[(i + 1) % 4]);
  }

  //look_down[i] = true if the frame is positioned below frame[i] and vice versa.
  bool look_down[4] = {}; //gotta zero all the values in look_down.

  for(int i = 0; i < 4; i++){
    //if the line is vertical then we don't need to compute anything
    if(!frame[i].vertical){
      //frame[i] holds the equation for the line thru shifted_frame[i] and shifted_frame[(i + 1)%4]
      //that means we have to check if one of the points not on the line frame[i] is below or above
      //said line.
      if(shifted_frame[(i + 2) % 4].y < (frame[i].m * shifted_frame[(i + 2) % 4].x + frame[i].b)){
        look_down[i] = 1;
      }
    }
  }

  R2Image::d_coordinates original_pix_pos;
  R2Pixel temp(0,0,0,1);
  double weight = 0;
  double total_weight = 0;

  bool in_frame = 1;
  for(int i = leftmost; i < rightmost; i++){
    for(int j = botmost; j < topmost; j++){
      in_frame = 1; //reset in_frame to be 1.
      for(int x = 0; x < 4; x++){
        //if the frame is vertical don't worry about it. 
        //you ain't got to touch it at all. the outer two for loops will deal with it.
        if(!frame[x].vertical){ 
          if(look_down[x]){ //if we should be looking down
            //this point should be below or on the line. if not then it's not in the frame.
            if(j > (frame[x].m * i + frame[x].b)){
              in_frame = 0;
            }
          }else{ //else we should be looking up
            //likewise, this point should be above or on the line, or else...
            if(j < (frame[x].m * i + frame[x].b)){
              in_frame = 0;
            } 
          }
        }
      }
      if(in_frame){
        //okay in here it's all about inverse warping.
        original_pix_pos = matrixMult(i, j, H);
        // if(original_pix_pos.x > width || original_pix_pos.x < 0 || original_pix_pos.y > height || original_pix_pos.y < 0){
        //   fprintf(stdout, "width and height of original: %d    %d \n", width, height);
        //   fprintf(stdout, "calculated position x and y: %d %d  \n", original_pix_pos.x, original_pix_pos.y);
        // }

        /*
        (floor, ceil)_____________(ceil, ceil)
              |                         |
              |             .(x,y)      |
              |                         |
        (floor, floor)____________(ceil, floor)


        */
        total_weight = 0;
        temp.Reset(0,0,0,0);

        //the weight for (ceil_x, ceil_y) is proportional to how far (x,y) is away from (floor, floor)
        weight = original_pix_pos.x - floor(original_pix_pos.x) + original_pix_pos.y - floor(original_pix_pos.y);
        temp = weight * this->Pixel(ceil(original_pix_pos.x), ceil(original_pix_pos.y));
        total_weight += weight;

        //the weight for (floor_x, floor_y) is proportional to how far (x,y) is away from (ceil, ceil)
        weight = ceil(original_pix_pos.x) - original_pix_pos.x + ceil(original_pix_pos.y) - original_pix_pos.y;
        temp += weight * this->Pixel(floor(original_pix_pos.x), floor(original_pix_pos.y));
        total_weight += weight;

        //the weight for (ceil_x, floor_y) is proportional to how far (x,y) is away from (floor, ceil)
        weight = original_pix_pos.x - floor(original_pix_pos.x) + ceil(original_pix_pos.y) - original_pix_pos.y;
        temp += weight * this->Pixel(ceil(original_pix_pos.x), floor(original_pix_pos.y));
        total_weight += weight;

        //the weight for (floor_x, ceil_y) is proportional to how far (x,y) is away from (ceil, floor)
        weight = ceil(original_pix_pos.x) - original_pix_pos.x + original_pix_pos.y - floor(original_pix_pos.y);
        temp += weight * this->Pixel(floor(original_pix_pos.x), ceil(original_pix_pos.y));
        total_weight += weight;

        //if the coordinates in original_pix_pos are whole numbers then we run the risk of dividing by zero
        if(total_weight != 0){ //and we never want to divide by zero.
          nextImage->Pixel(i, j) = temp/total_weight; 
        }else{
          nextImage->Pixel(i, j) = this->Pixel(original_pix_pos.x, original_pix_pos.y);
        }
      }
    }
  }

}





////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
  // Read pixel values
  int red, green, blue;
  if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
    fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
    fclose(fp);
    return 0;
  }

  // Assign pixel values
  double r = (double) red / max_value;
  double g = (double) green / max_value;
  double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


  

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width;  /* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;   /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;   /* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
  
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






