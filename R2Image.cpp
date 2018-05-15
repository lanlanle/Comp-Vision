// Source file for image class



// Include files 
#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <math.h> 
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include "algorithm"
#include <iomanip>
#include <limits>
using namespace std;



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
	// R2Point p1(1.2,3.5);
	// R2Point p2(2.1,2.2);
	// R2Point p3(0.2,1.6);
	// R2Point p4(0.0,0.5);
	// R2Point p5(-0.2,4.2);

	// // build the 5x6 matrix of equations
	// double** linEquations = dmatrix(1,5,1,6);

	// linEquations[1][1] = p1[0]*p1[0];
	// linEquations[1][2] = p1[0]*p1[1];
	// linEquations[1][3] = p1[1]*p1[1];
	// linEquations[1][4] = p1[0];
	// linEquations[1][5] = p1[1];
	// linEquations[1][6] = 1.0;

	// linEquations[2][1] = p2[0]*p2[0];
	// linEquations[2][2] = p2[0]*p2[1];
	// linEquations[2][3] = p2[1]*p2[1];
	// linEquations[2][4] = p2[0];
	// linEquations[2][5] = p2[1];
	// linEquations[2][6] = 1.0;

	// linEquations[3][1] = p3[0]*p3[0];
	// linEquations[3][2] = p3[0]*p3[1];
	// linEquations[3][3] = p3[1]*p3[1];
	// linEquations[3][4] = p3[0];
	// linEquations[3][5] = p3[1];
	// linEquations[3][6] = 1.0;
	
	// linEquations[4][1] = p4[0]*p4[0];
	// linEquations[4][2] = p4[0]*p4[1];
	// linEquations[4][3] = p4[1]*p4[1];
	// linEquations[4][4] = p4[0];
	// linEquations[4][5] = p4[1];
	// linEquations[4][6] = 1.0;

	// linEquations[5][1] = p5[0]*p5[0];
	// linEquations[5][2] = p5[0]*p5[1];
	// linEquations[5][3] = p5[1]*p5[1];
	// linEquations[5][4] = p5[0];
	// linEquations[5][5] = p5[1];
	// linEquations[5][6] = 1.0;

	// printf("\n Fitting a conic to five points:\n");
	// printf("Point #1: %f,%f\n",p1[0],p1[1]);
	// printf("Point #2: %f,%f\n",p2[0],p2[1]);
	// printf("Point #3: %f,%f\n",p3[0],p3[1]);
	// printf("Point #4: %f,%f\n",p4[0],p4[1]);
	// printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// // compute the SVD
	// double singularValues[7]; // 1..6
	// double** nullspaceMatrix = dmatrix(1,6,1,6);
	// svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// // get the result
	// printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// // find the smallest singular value:
	// int smallestIndex = 1;
	// for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	// printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// // make sure the solution is correct:
	// printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
	// 									p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
	// 									p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
	// 									p1[0]*nullspaceMatrix[4][smallestIndex] + 
	// 									p1[1]*nullspaceMatrix[5][smallestIndex] + 
	// 									nullspaceMatrix[6][smallestIndex]);

	// printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
	// 									p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
	// 									p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
	// 									p2[0]*nullspaceMatrix[4][smallestIndex] + 
	// 									p2[1]*nullspaceMatrix[5][smallestIndex] + 
	// 									nullspaceMatrix[6][smallestIndex]);

	// printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
	// 									p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
	// 									p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
	// 									p3[0]*nullspaceMatrix[4][smallestIndex] + 
	// 									p3[1]*nullspaceMatrix[5][smallestIndex] + 
	// 									nullspaceMatrix[6][smallestIndex]);

	// printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
	// 									p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
	// 									p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
	// 									p4[0]*nullspaceMatrix[4][smallestIndex] + 
	// 									p4[1]*nullspaceMatrix[5][smallestIndex] + 
	// 									nullspaceMatrix[6][smallestIndex]);

	// printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
	// 									p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
	// 									p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
	// 									p5[0]*nullspaceMatrix[4][smallestIndex] + 
	// 									p5[1]*nullspaceMatrix[5][smallestIndex] + 
	// 									nullspaceMatrix[6][smallestIndex]);

	// R2Point test_point(0.34,-2.8);

	// printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
	// 										test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
	// 										test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
	// 										test_point[0]*nullspaceMatrix[4][smallestIndex] + 
	// 										test_point[1]*nullspaceMatrix[5][smallestIndex] + 
	// 										nullspaceMatrix[6][smallestIndex]);

	// return;	

 


  R2Point a1(0,0);
  R2Point a2(1,0);
  R2Point a3(1,1);
  R2Point a4(0,1);

  vector<R2Point> original_vector;

  original_vector.push_back(a1);
  original_vector.push_back(a2);
  original_vector.push_back(a3);
  original_vector.push_back(a4);

  // R2Point b1(1,0);
  // R2Point b2(2,0);
  // R2Point b3(2,1);
  // R2Point b4(1,1);

  R2Point b1(1,2);
  R2Point b2(1,1);
  R2Point b3(3,1);
  R2Point b4(3,2);

  vector<R2Point> transformed_vector;

  transformed_vector.push_back(b1);
  transformed_vector.push_back(b2);
  transformed_vector.push_back(b3);
  transformed_vector.push_back(b4);


  // construct the A matrix 
  // vector < vector <double> > A;


  double** linEquations = dmatrix(1,8,1,9);

  vector <vector <double> > A;

  for (int i = 0; i<original_vector.size();i++) {
    vector <double> equation1;

    linEquations[i*2+1][1] = 0.0;
    linEquations[i*2+1][2] = 0.0;
    linEquations[i*2+1][3] = 0.0;
    linEquations[i*2+1][4] = -original_vector[i][0];
    linEquations[i*2+1][5] = -original_vector[i][1];
    linEquations[i*2+1][6] = -1.0;
    linEquations[i*2+1][7] = transformed_vector[i][1]*original_vector[i][0];
    linEquations[i*2+1][8] = transformed_vector[i][1]*original_vector[i][1];
    linEquations[i*2+1][9] = transformed_vector[i][1];

      equation1.push_back(0.0);
      equation1.push_back(0.0);
      equation1.push_back(0.0);
      equation1.push_back(-original_vector[i][0]);
      equation1.push_back(-original_vector[i][1]);
      equation1.push_back(-1.0);
      equation1.push_back(transformed_vector[i][1]*original_vector[i][0]);
      equation1.push_back(transformed_vector[i][1]*original_vector[i][0]);
      equation1.push_back(transformed_vector[i][1]);

      A.push_back(equation1);


    vector<double> equation2;

    linEquations[i*2+2][1] = original_vector[i][0];
    linEquations[i*2+2][2] = original_vector[i][1];
    linEquations[i*2+2][3] = 1.0;
    linEquations[i*2+2][4] = 0.0;
    linEquations[i*2+2][5] = 0.0;
    linEquations[i*2+2][6] = 0.0;
    linEquations[i*2+2][7] = -transformed_vector[i][0]*original_vector[i][0];
    linEquations[i*2+2][8] = -transformed_vector[i][0]*original_vector[i][1];
    linEquations[i*2+2][9] = -transformed_vector[i][0];

    equation2.push_back(original_vector[i][0]);
    equation2.push_back( original_vector[i][1]);
    equation2.push_back(1.0);
    equation2.push_back(0.0);
    equation2.push_back(0.0);
    equation2.push_back(0.0);
    equation2.push_back( -transformed_vector[i][0]*original_vector[i][0]);
    equation2.push_back( -transformed_vector[i][0]*original_vector[i][1]);
    equation2.push_back( -transformed_vector[i][0]);
    A.push_back(equation2);



  }

  //print A matrix

  cout << "A MATRIX: " << endl;

  for (int i = 0; i < 8; i++){
    cout<< A[i][0]<< " "<< A[i][1]<< " "<< A[i][2]<< " "<< A[i][3]<< " "<< A[i][4]<< " "<< A[i][5]<< " "<< A[i][6]<< " "<< A[i][7]<< " "<< A[i][8]<< endl;
  }
  cout << " " << endl;


  // // compute the SVD
  double singularValues[10]; // 1..6
  double** nullspaceMatrix = dmatrix(1,9,1,9);
  svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);


  for (int i = 1; i<=9; i++) {
    cout << "Singular Values " << i << ": " << setprecision(5) << singularValues[i] << endl;
  }
   cout << " " << endl;

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<=9;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  // printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

   for (int i = 1; i<=9; i++) {
    cout << "The H value " << i << ": " << setprecision(5) << nullspaceMatrix[i][smallestIndex] << endl;
  }
   cout << " " << endl;


  //make sure the algs are correct 
  for (int i = 0; i < 8; i++) {
    double result = A[i][0]*nullspaceMatrix[1][smallestIndex]+
                    A[i][1]*nullspaceMatrix[2][smallestIndex]+
                    A[i][2]*nullspaceMatrix[3][smallestIndex]+
                    A[i][3]*nullspaceMatrix[4][smallestIndex]+
                    A[i][4]*nullspaceMatrix[5][smallestIndex]+
                    A[i][5]*nullspaceMatrix[6][smallestIndex]+
                    A[i][6]*nullspaceMatrix[7][smallestIndex]+
                    A[i][7]*nullspaceMatrix[8][smallestIndex]+
                    A[i][8]*nullspaceMatrix[9][smallestIndex];

      cout << "Equation " << setprecision(5) << i+1 << ": " << result << endl;

  }



}


void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if(x0>x1)
  {
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
   if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
           // note that this division needs to be done in a way that preserves the fractional part
     int y = y0;
     for(int x=x0;x<=x1;x++)
   {
     Pixel(x,y).Reset(r,g,b,0.2);
         error = error + deltaerr;
         if(error>=0.5)
     {
       if(deltay>0) y = y + 1;
       else y = y - 1;

             error = error - 1.0;
     }
   }
   if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
   {
     for(int x=x0-3;x<=x0+3;x++)
     {
       for(int y=y0-3;y<=y0+3;y++)
       {
         Pixel(x,y).Reset(r,g,b,0.2);
       }
     }
   }
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
SobelX(void)
{

  //change to grayscale
  // for (int i = 0; i < width; i++) {
  //   for (int j = 0;  j < height; j++) {
  //     double r =Pixel(i,j).Red();
  //     double g =Pixel(i,j).Green();
  //     double b =Pixel(i,j).Blue();
  //     double a = Pixel(i,j).Alpha();
  //     double grayscale =0.2989*r + 0.5870*g + 0.1140*b; 

  //     // if(grayscale>1) grayscale = 1;
  //     // if(grayscale< 0) grayscale = 0;

  //     Pixel(i,j).Reset(grayscale,grayscale,grayscale,a);
  //   }
  // }
  //make a copy of the image 
  R2Image temp_img(*this);
	// Apply the Sobel oprator to the image in X direction
  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      
      double newRed = temp_img.Pixel(i-1,j-1).Red()*(-1)+temp_img.Pixel(i-1,j).Red()*(-2)+temp_img.Pixel(i-1,j+1).Red()*(-1)+temp_img.Pixel(i+1,j-1).Red()+temp_img.Pixel(i+1,j).Red()*2+temp_img.Pixel(i+1,j+1).Red();
      double newGreen = temp_img.Pixel(i-1,j-1).Green()*(-1)+temp_img.Pixel(i-1,j).Green()*(-2)+temp_img.Pixel(i-1,j+1).Green()*(-1)+temp_img.Pixel(i+1,j-1).Green()+temp_img.Pixel(i+1,j).Green()*2+temp_img.Pixel(i+1,j+1).Green();
      double newBlue = temp_img.Pixel(i-1,j-1).Blue()*(-1)+temp_img.Pixel(i-1,j).Blue()*(-2)+temp_img.Pixel(i-1,j+1).Blue()*(-1)+temp_img.Pixel(i+1,j-1).Blue()+temp_img.Pixel(i+1,j).Blue()*2+temp_img.Pixel(i+1,j+1).Blue();
      double alpha = Pixel(i,j).Alpha();
      Pixel(i,j).Reset(newRed,newGreen,newBlue,alpha);
      // Pixel(i,j).Clamp(1);

    }
  }
}

void R2Image::
SobelY(void)
{
  //change img to grayscale
  // for (int i = 0; i < width; i++) {
  //   for (int j = 0;  j < height; j++) {
  //     double r =Pixel(i,j).Red();
  //     double g =Pixel(i,j).Green();
  //     double b =Pixel(i,j).Blue();
  //     double a = Pixel(i,j).Alpha();
  //     double grayscale =0.2989*r + 0.5870*g + 0.1140*b; 

  //     // if(grayscale>1) grayscale = 1;
  //     // if(grayscale< 0) grayscale = 0;

  //     Pixel(i,j).Reset(grayscale,grayscale,grayscale,a);
  //     // Pixel(i,j).Clamp(1);
  //   }
  // }
	//how to make a copy of the image 
  R2Image temp_img(*this);
  // Apply the Sobel oprator to the image in  direction
  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      double newRed = temp_img.Pixel(i-1,j-1).Red()*(-1)+temp_img.Pixel(i,j-1).Red()*(-2)+temp_img.Pixel(i+1,j-1).Red()*(-1)+temp_img.Pixel(i-1,j+1).Red()+temp_img.Pixel(i,j+1).Red()*2+temp_img.Pixel(i+1,j+1).Red();
      double newGreen = temp_img.Pixel(i-1,j-1).Green()*(-1)+temp_img.Pixel(i,j-1).Green()*(-2)+temp_img.Pixel(i+1,j-1).Green()*(-1)+temp_img.Pixel(i-1,j+1).Green()+temp_img.Pixel(i,j+1).Green()*2+temp_img.Pixel(i+1,j+1).Green();
      double newBlue = temp_img.Pixel(i-1,j-1).Blue()*(-1)+temp_img.Pixel(i,j-1).Blue()*(-2)+temp_img.Pixel(i+1,j-1).Blue()*(-1)+temp_img.Pixel(i-1,j+1).Blue()+temp_img.Pixel(i,j+1).Blue()*2+temp_img.Pixel(i+1,j+1).Blue();
      double alpha = Pixel(i,j).Alpha();

      Pixel(i,j).Reset(newRed,newGreen,newBlue,alpha);
      // Pixel(i,j).Clamp(1);

    }
  }
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}


void R2Image::
Median(void)
{
  // Apply the Median filter to the image

  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  // fprintf(stderr, "Median not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      double r =Pixel(i,j).Red();
      double g =Pixel(i,j).Green();
      double b =Pixel(i,j).Blue();
      double a = Pixel(i,j).Alpha();
      double grayscale =0.2989*r + 0.5870*g + 0.1140*b; 
      double newRed = r*factor-grayscale*(factor-1);
      double newGreen = g*factor-grayscale*(factor-1);
      double newBlue = b*factor-grayscale*(factor-1);

      if(newRed>1) newRed = 1;
      if(newGreen >1 ) newGreen = 1;
      if(newBlue >1) newBlue = 1;

      if(newRed< 0) newRed = 0;
      if(newGreen < 0) newGreen = 0;
      if(newBlue < 0) newBlue = 0;
      Pixel(i,j).Reset(newRed,newGreen,newBlue,a);
    }
  }
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{

  //reset 
  R2Image temp_img(*this);
  // Gaussian blur of the image. Separable solution is preferred

  
  int size = 6*sigma+1;
  int half_size = size/2;
  vector<double> kernel;
  //normalize 
  // double sum = 0.0;


  // set up the kernel

  for (int i = 0; i < size; i++){
    int x = i-half_size;
    double value =(exp(-(x*x)/(2*sigma*sigma)))/sqrt(2*3.14159265*sigma*sigma);
    kernel.push_back(value);
    // sum+=value;
  }


  //filter horizontally 
  for (int rh = 0; rh < width; rh++) {
    for (int ch = 0; ch < height; ch++) {
      //new value 
      R2Pixel* temp_pixel = new R2Pixel();
      for (int pos = -(half_size); pos <= half_size; pos++) {
          int col_pos= ch+pos;
          if (col_pos<0) col_pos=0;
          if (col_pos>=height)col_pos=height-1;          
          *temp_pixel += kernel[pos+half_size]*temp_img.Pixel(rh,col_pos);
          
      }
      Pixel(rh,ch)= *temp_pixel;
      // Pixel(rh,ch).Clamp(1);
      
    }
  }
  R2Image temp_horizontal(*this);

  //filter vertically 
  for (int cv = 0; cv < height;cv++ ){
    for (int rv= 0; rv < width; rv++){
      R2Pixel* temp_pixel = new R2Pixel();
      for (int pos = -(half_size); pos <= half_size; pos++) {
          int row_pos= rv+pos;
          if (row_pos<0) row_pos=0;
          if (row_pos>= width)row_pos = width-1;
          *temp_pixel += kernel[pos+half_size]*temp_horizontal.Pixel(row_pos,cv);
          
      }
      Pixel(rv,cv)= *temp_pixel;
      // Pixel(rv,cv).Clamp(1);
    }
  }
}



vector<vector<double> > R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

  R2Image original(*this);

  R2Image Img1(*this);
  Img1.SobelX();
  R2Image Img2(*this);
  Img2.SobelY();
  R2Image Img3(*this);
  Img3.SobelX();


  for (int i = 0; i < width; i++) {
      for (int j = 0 ; j< height; j++) {
          
        Img3.Pixel(i,j) = Img3.Pixel(i,j)* Img2.Pixel(i,j);
        Img1.Pixel(i,j) = Img1.Pixel(i,j)* Img1.Pixel(i,j);
        Img2.Pixel(i,j) = Img2.Pixel(i,j)* Img2.Pixel(i,j);

      }
  }


  // set up the image with 

  Img1.Blur(sigma);
  Img2.Blur(sigma);
  Img3.Blur(sigma);

  for (int i=0; i < width; i++){
    for (int j=0; j<height;j++) {
      R2Pixel* temp_pixel = new R2Pixel();
      *temp_pixel = Img1.Pixel(i,j)*Img2.Pixel(i,j)-Img3.Pixel(i,j)*Img3.Pixel(i,j)-0.04*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));

      double r = temp_pixel->Red()+0.5;
      double g = temp_pixel->Green()+0.5;
      double b = temp_pixel->Blue()+0.5;
      double a = temp_pixel->Alpha();

      R2Pixel* adjusted_pixel = new R2Pixel(r,g,b,a);
      Pixel(i,j) = *adjusted_pixel;
      Pixel(i,j).Clamp();
      
    }
  }

  //feature description

  vector< vector<double> > corners;  

  //calculate Harris Score for every pixel
  for (int i=5; i<width-5; i++){
    for (int j=5; j<height-5;j++){
      vector<double> item;

      double harris_score = Pixel(i,j).Red()+Pixel(i,j).Green()+Pixel(i,j).Blue();
      item.push_back(harris_score);
      item.push_back(i);
      item.push_back(j);
      corners.push_back(item);

    }
  }
  //Sort to get 150 pixels with highest Harris score 

  sort(corners.begin(),corners.end(),[](const vector<double>& a, const std::vector<double>& b) {
  return a[0] > b[0];
});

  // choose 150 separate features 
  vector< vector<double> > displayed_pixels;
  int counter= 0; // how many items already in the displayed_pixel list
  int currentItem = 0;// items being looked at 

  while (counter<=150 && currentItem < corners.size()) {
      bool isOutOfBound=true;
      for (int i=0; i<displayed_pixels.size(); i++) {
        // if the current pixel being looked at is not out of bound 
        if (corners[currentItem][1] <= displayed_pixels[i][1]+10 && 
            corners[currentItem][1] >= displayed_pixels[i][1]-10 &&
            corners[currentItem][2] <= displayed_pixels[i][2]+10 && 
            corners[currentItem][2] >= displayed_pixels[i][2]-10) {
            isOutOfBound = false;
            break;
        }
      }
     
      if (isOutOfBound){
        displayed_pixels.push_back(corners[currentItem]);
        counter++;
      } 
      currentItem++;
  }

  // //draw the pixel 

  // for (int pixel = 0; pixel< displayed_pixels.size();pixel++) {
  //   //draw scquares around that feature 
  //   for (int i = -5; i <=5; i++){
  //     original.Pixel(displayed_pixels[pixel][1]+i,displayed_pixels[pixel][2]-5).Reset(1,0,0,1);
  //     original.Pixel(displayed_pixels[pixel][1]+i,displayed_pixels[pixel][2]+5).Reset(1,0,0,1);
  //     original.Pixel(displayed_pixels[pixel][1]-5,displayed_pixels[pixel][2]+i).Reset(1,0,0,1);
  //     original.Pixel(displayed_pixels[pixel][1]+5,displayed_pixels[pixel][2]+i).Reset(1,0,0,1);
  //   }
    
  // }
  // Clamp the pixels

  for (int i =0; i < width; i++){
    for (int j = 0; j< height; j++) {
      Pixel(i,j)= original.Pixel(i,j);
      Pixel(i,j).Clamp();
    }
  }


  return displayed_pixels;
}




void R2Image::
Sharpen()
{
  
  R2Image temp_img(*this);
 // Sharpen an image using a linear filter. Use a kernel of your choosing.
  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      
      double newRed = temp_img.Pixel(i,j-1).Red()*(-1)+temp_img.Pixel(i-1,j).Red()*(-1)+temp_img.Pixel(i,j+1).Red()*(-1)+temp_img.Pixel(i+1,j).Red()*(-1)+temp_img.Pixel(i,j).Red()*5;
      double newGreen = temp_img.Pixel(i,j-1).Green()*(-1)+temp_img.Pixel(i-1,j).Green()*(-1)+temp_img.Pixel(i,j+1).Green()*(-1)+temp_img.Pixel(i+1,j).Green()*(-1)+temp_img.Pixel(i,j).Green()*5;
      double newBlue = temp_img.Pixel(i,j-1).Blue()*(-1)+temp_img.Pixel(i-1,j).Blue()*(-1)+temp_img.Pixel(i,j+1).Blue()*(-1)+temp_img.Pixel(i+1,j).Blue()*(-1)+temp_img.Pixel(i,j).Blue()*5;
      double alpha = Pixel(i,j).Alpha();

      Pixel(i,j).Reset(newRed,newGreen,newBlue,alpha);
      Pixel(i,j).Clamp(1);

    }
  }
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
	// into this image with a 50% opacity.

  vector<vector<double>> imgA_features = Harris(2.0);

  vector<vector<double> > matching_features;

  
  int searchAreaHalfWidth = otherImage->width * 0.01;
  int searchAreaHalfHeight = otherImage->height * 0.01;
   cout << searchAreaHalfWidth << endl;
   cout << searchAreaHalfHeight << endl;

  //search for all feature in vector
  for (int feature = 0; feature < imgA_features.size(); feature++) {
    double minSum = 100000; 
    double minX = 0;
    double minY = 0;
    for (int i = -searchAreaHalfWidth; i <=searchAreaHalfWidth;i++) {
      for (int j = -searchAreaHalfHeight; j <= searchAreaHalfHeight; j++) {

        if (imgA_features[feature][1]+i >= 0 && imgA_features[feature][1]+i< width &&
                imgA_features[feature][2]+j>=0 && imgA_features[feature][2]+j<height){

                 // start a window
                  R2Pixel *ssd = new R2Pixel();
                  for (int k = -3; k <=3; k++) {
                    for (int l = -3; l <=3; l++) {
                      if (imgA_features[feature][1]+k >= 0 && imgA_features[feature][1]+k < otherImage->width &&
                          imgA_features[feature][2]+l >=0 && imgA_features[feature][2]+l < otherImage->height &&
                          imgA_features[feature][1]+i+k >= 0 && imgA_features[feature][1]+i+k < otherImage->width &&
                          imgA_features[feature][2]+j+l>=0 && imgA_features[feature][2]+j+l<otherImage->height) 
             
                          *ssd += (Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-otherImage->Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l))*(Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-otherImage->Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l));          
                      
                    }
                  }
                  double sum = ssd->Red()+ssd->Green()+ssd->Blue();

                  if (sum < minSum) {
                    minSum = sum; 
                    minX = imgA_features[feature][1]+i;
                    minY = imgA_features[feature][2]+j;
                  
                  }

        }
            
      }
    }

    vector<double> item;
    item.push_back(minX);
    item.push_back(minY);
    matching_features.push_back(item);
  }

  //draw the matching features 
  for (int i = 0; i< matching_features.size(); i++){
    line(matching_features[i][0],matching_features[i][0],matching_features[i][1],matching_features[i][1],0,1,0);
  }



  // use random to find a vector item 
  // int max_inlier = 0;
  // vector < vector<double> > max_inlier_set;
  // vector < vector<double> > min_outlier_set; 

  // for (int trial=0; trial < 200; trial++) {
  //   int randomIndex = rand() % matching_features.size();

  
   
  //   //values for the  sample vector 
  //   double sampleX = imgA_features[randomIndex][1]-matching_features[randomIndex][0];
  //   double sampleY = imgA_features[randomIndex][2]-matching_features[randomIndex][1];
    

  //   vector<vector<double> > inlier;
  //   vector <vector<double> > outlier;        
  //   for (int v = 0; v < matching_features.size();v++){

  //       //values for the current vector 
  //       double currentX = imgA_features[v][1]-matching_features[v][0];
  //       double currentY = imgA_features[v][2]-matching_features[v][1];

  //       double difference = sqrt(pow(currentX-sampleX,2)+pow(currentY-sampleY,2));
      

        


  //       vector<double> motion_vector;
  //       motion_vector.push_back(imgA_features[v][1]);
  //       motion_vector.push_back(matching_features[v][0]);
  //       motion_vector.push_back(imgA_features[v][2]);
  //       motion_vector.push_back(matching_features[v][1]);

  //       if (difference< 4.0) inlier.push_back(motion_vector);
  //       else outlier.push_back(motion_vector);

  //   }
  //   if (inlier.size() > max_inlier){
  //       max_inlier = inlier.size();
  //       max_inlier_set.clear();
  //       min_outlier_set.clear();

  //       for (int a= 0; a < inlier.size(); a++) max_inlier_set.push_back(inlier[a]);
  //       for (int b= 0; b < outlier.size(); b++) min_outlier_set.push_back(outlier[b]);

  //   }

  // }
  // cout <<"INLIER: " << max_inlier << endl;

  // //draw inlier set

  // for (int a = 0 ; a < max_inlier_set.size(); a++) {
    
  //       otherImage->line(max_inlier_set[a][0],max_inlier_set[a][1],max_inlier_set[a][2],max_inlier_set[a][3],0,1,0);
  // }

  // // draw outlier 
  // for (int a = 0 ; a < min_outlier_set.size(); a++) {
    
  //       otherImage->line(min_outlier_set[a][0],min_outlier_set[a][1],min_outlier_set[a][2],min_outlier_set[a][3],1,0,0);
  // }
}

vector<vector<double>> imgA_features; 


void R2Image::FirstFrameProcessing(){
    localFirstImage = new R2Image(*this);
    imgA_features= localFirstImage->Harris(2.0);
    TrackFeatures(this);
}



void R2Image::FrameProcessing(R2Image *otherImage){
  otherImage->TrackFeatures(this);
}

    


void R2Image::TrackFeatures(R2Image *otherImage){
  


  vector<vector<double> > matching_features;


  
  int searchAreaHalfWidth = width * 0.1;
  int searchAreaHalfHeight = height * 0.1;
   cout << searchAreaHalfWidth << endl;
   cout << searchAreaHalfHeight << endl;

  //search for all feature in vector
  for (int feature = 0; feature < imgA_features.size(); feature++) {
    double minSum = 100000; 
    double minX = 0;
    double minY = 0;
    for (int i = -searchAreaHalfWidth; i <=searchAreaHalfWidth;i++) {
      for (int j = -searchAreaHalfHeight; j <= searchAreaHalfHeight; j++) {

        if (imgA_features[feature][1]+i >= 0 && imgA_features[feature][1]+i<width &&
                imgA_features[feature][2]+j>=0 && imgA_features[feature][2]+j<height){

                 // start a window
                  R2Pixel *ssd = new R2Pixel();
                  for (int k = -3; k <=3; k++) {
                    for (int l = -3; l <=3; l++) {
                      if (imgA_features[feature][1]+k >= 0 && imgA_features[feature][1]+k < width &&
                          imgA_features[feature][2]+l >=0 && imgA_features[feature][2]+l < height &&
                          imgA_features[feature][1]+i+k >= 0 && imgA_features[feature][1]+i+k < width &&
                          imgA_features[feature][2]+j+l>=0 && imgA_features[feature][2]+j+l< height) 
             
                         *ssd += (Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-otherImage->Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l))*(Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-otherImage->Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l));          
                      
                      
                    }
                  }
                  double sum = ssd->Red()+ssd->Green()+ssd->Blue();

                  if (sum < minSum) {
                    minSum = sum; 
                    minX = imgA_features[feature][1]-i;
                    minY = imgA_features[feature][2]-j;
                  
                  }

        }
            
      }
    }
       cout << "difference: "<< minSum << endl;
    
      vector<double> item;
      item.push_back(minX);
      item.push_back(minY);
      matching_features.push_back(item);
    
  }

  //draw the matching features 
  for (int i = 0; i< matching_features.size(); i++){
   line(matching_features[i][0],matching_features[i][0],matching_features[i][1],matching_features[i][1],0,1,0);
  }

  imgA_features.clear();

  for(int i = 0; i < matching_features.size();i++){
    vector<double> object;
    object.push_back(0.0);
    // cout << "matching features: "<< matching_features[i][0]<< " "<<matching_features[i][1]<<endl;
    object.push_back(matching_features[i][0]);
    object.push_back(matching_features[i][1]);

    imgA_features.push_back(object);
  }

   

}


vector<double> R2Image::svdHelper(vector<vector <double> > A, vector<vector <double> > B, int run){

   // instead of one, choose 4 here 
 
  vector<R2Point> original_vector;
  vector<R2Point> transformed_vector;


  if (run==4){

    vector <double> random_vector;
    for (int i = 0 ; i< run; i++) {
      
      int randomIndex = rand() % B.size();

      random_vector.push_back(randomIndex);
      // cout << randomIndex << endl;
    }
    for (int i = 0; i < random_vector.size(); i++) {
      R2Point a(A[random_vector[i]][1],A[random_vector[i]][2]);
      original_vector.push_back(a);
    }
    

      for (int i = 0; i < random_vector.size(); i++) {
        R2Point b(B[random_vector[i]][0],B[random_vector[i]][1]);
        transformed_vector.push_back(b);
      }

  } else {
      for (int i = 0 ; i<A.size(); i++){
        R2Point a(A[i][1],A[i][2]);
        original_vector.push_back(a);
        R2Point b(B[i][0],B[i][1]);
        transformed_vector.push_back(b);
      } 
  }

    
    


    //all the inlier points 
    double** linEquations = dmatrix(1,run*2,1,9);

    for (int i = 0; i<original_vector.size();i++) {

      linEquations[i*2+1][1] = 0.0;
      linEquations[i*2+1][2] = 0.0;
      linEquations[i*2+1][3] = 0.0;
      linEquations[i*2+1][4] = -original_vector[i][0];
      linEquations[i*2+1][5] = -original_vector[i][1];
      linEquations[i*2+1][6] = -1.0;
      linEquations[i*2+1][7] = transformed_vector[i][1]*original_vector[i][0];
      linEquations[i*2+1][8] = transformed_vector[i][1]*original_vector[i][1];
      linEquations[i*2+1][9] = transformed_vector[i][1];

      linEquations[i*2+2][1] = original_vector[i][0];
      linEquations[i*2+2][2] = original_vector[i][1];
      linEquations[i*2+2][3] = 1.0;
      linEquations[i*2+2][4] = 0.0;
      linEquations[i*2+2][5] = 0.0;
      linEquations[i*2+2][6] = 0.0;
      linEquations[i*2+2][7] = -transformed_vector[i][0]*original_vector[i][0];
      linEquations[i*2+2][8] = -transformed_vector[i][0]*original_vector[i][1];
      linEquations[i*2+2][9] = -transformed_vector[i][0];

    }




  // // compute the SVD
  double singularValues[10]; // 1..6
  double** nullspaceMatrix = dmatrix(1,9,1,9);
  svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);


  

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<=9;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;


  //WE GET H VALUES HERE! 
    vector<double> H_values;

  for (int i = 1; i<=9; i++) {
        // cout << "The H value " << i << ": " << nullspaceMatrix[i][smallestIndex] << endl;
        H_values.push_back(nullspaceMatrix[i][smallestIndex]);
  }
  return H_values;
   

}


  


void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
	// fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	// return;



  //PRODUCE 150 VECTOR TRANSLATION 
   R2Image temp(*this);

   R2Image newA(*otherImage);


  vector<vector<double>> imgA_features = otherImage->Harris(2.0);

  vector<vector<double> > matching_features;

  
  int searchAreaHalfWidth = width * 0.05;
  int searchAreaHalfHeight = height * 0.05;
   cout << searchAreaHalfWidth << endl;
   cout << searchAreaHalfHeight << endl;

  //search for all feature in vector
  for (int feature = 0; feature < imgA_features.size(); feature++) {
    double minSum = 100000; 
    double minX = 0;
    double minY = 0;
    for (int i = -searchAreaHalfWidth; i <= searchAreaHalfWidth;i++) {
      for (int j = -searchAreaHalfHeight; j <= searchAreaHalfHeight; j++) {

        if (imgA_features[feature][1]+i >= 0 && imgA_features[feature][1]+i< width &&
                imgA_features[feature][2]+j>=0 && imgA_features[feature][2]+j<height){

                 // start a window
                  R2Pixel *ssd = new R2Pixel();
                  for (int k = -3; k <=3; k++) {
                    for (int l = -3; l <=3; l++) {
                      if (imgA_features[feature][1]+k >= 0 && imgA_features[feature][1]+k < width &&
                          imgA_features[feature][2]+l >=0 && imgA_features[feature][2]+l < height &&
                          imgA_features[feature][1]+i+k >= 0 && imgA_features[feature][1]+i+k < width &&
                          imgA_features[feature][2]+j+l>=0 && imgA_features[feature][2]+j+l<height) 
             
                          *ssd += (otherImage->Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l))*(otherImage->Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l));          
                      
                    }
                  }
                  double sum = ssd->Red()+ssd->Green()+ssd->Blue();

                  if (sum < minSum) {
                    minSum = sum; 
                    minX = imgA_features[feature][1]+i;
                    minY = imgA_features[feature][2]+j;
                  
                  }
        }  
        
      }
    }

    vector<double> item;
    item.push_back(minX);
    item.push_back(minY);
    matching_features.push_back(item);
  }



  // use random to find a vector item 
  int max_inlier = 0;
  vector < vector<double> > max_inlier_orig_set;
  vector < vector<double> > max_inlier_trans_set; 
  vector <double> temp_H_matrix; 




  // create the model of H matrix first here 



  for (int trial=0; trial < 1000; trial++) {

  vector<double> hValuesBA = svdHelper(imgA_features,matching_features,4);
  cout << " " << endl;


  //COMPARE TO SEE IF THE TRANSLATION IS CORRECT 

  vector<vector<double> > inlier_orig;
  vector <vector<double> > inlier_trans;   
  for (int i = 0; i < matching_features.size(); i++)   {
        
        double xB = hValuesBA[0]*imgA_features[i][1]+hValuesBA[1]*imgA_features[i][2]+hValuesBA[2];
        double yB = hValuesBA[3]*imgA_features[i][1]+hValuesBA[4]*imgA_features[i][2]+hValuesBA[5];
        double wB = hValuesBA[6]*imgA_features[i][1]+hValuesBA[7]*imgA_features[i][2]+hValuesBA[8];
       
        
        double scaled_xB = xB/wB;
        double scaled_yB = yB/wB; 

       
        double distance = sqrt(pow(matching_features[i][0]-scaled_xB,2)+pow(matching_features[i][1]-scaled_yB,2));


        vector<double> motion_vector1;
        motion_vector1.push_back(imgA_features[i][0]);
        motion_vector1.push_back(imgA_features[i][1]);
        motion_vector1.push_back(imgA_features[i][2]);
        

        vector<double> motion_vector2;
        motion_vector2.push_back(matching_features[i][0]);
        motion_vector2.push_back(matching_features[i][1]);


        if ( distance <=5) {
          inlier_orig.push_back(motion_vector1);
          inlier_trans.push_back(motion_vector2);
        }
        
    }
    
    if (inlier_orig.size() > max_inlier){
        max_inlier = inlier_orig.size();
        max_inlier_orig_set.clear();
        max_inlier_trans_set.clear();
        temp_H_matrix.clear();
    
        for (int a= 0; a < inlier_orig.size(); a++){
           max_inlier_orig_set.push_back(inlier_orig[a]);
           max_inlier_trans_set.push_back(inlier_trans[a]);
        }  
        for (int d = 0; d<9; d++) temp_H_matrix.push_back(hValuesBA[d]);

    }

  }

  cout << " MATRIX_VALUES" << endl;
  for (int i = 0; i < 9; i++) cout << temp_H_matrix[i] << " ";
  cout << " "<< endl;

  cout <<"INLIER: " << max_inlier << endl;
  cout << "REAL SIZE" << max_inlier_orig_set.size() << endl;

  vector<vector<double> > vecA;
  vector<vector<double> > vecB;
  for (int i = 0; i < max_inlier; i++){
    vecA.push_back(max_inlier_orig_set[i]);
    vecB.push_back(max_inlier_trans_set[i]);
  }

  vector<double> best_H_matrix= svdHelper(max_inlier_orig_set,max_inlier_trans_set,max_inlier_orig_set.size());


  //use best inverse matrix to reconstruct a new image 
  cout << " MATRIX_VALUES" << endl;
  for (int i = 0; i < 9; i++) cout << best_H_matrix[i] << " ";
  cout << " "<< endl;


  //compute matrix inverse
  vector<vector <double>> mat; 
  for (int i= 0; i<3;i++) {
    vector<double> line;
    for (int j= 0;j<3; j++) {
      line.push_back(best_H_matrix[3*i+j]);
    }
    mat.push_back(line);
  }


  vector<double>inverse_mat;
  float determinant = mat[0][0]*mat[1][1]*mat[2][2]+mat[0][1]*mat[1][2]*mat[2][0]+mat[0][2]*mat[1][0]*mat[2][1]
                      -(mat[0][2]*mat[1][1]*mat[2][0]+mat[0][0]*mat[1][2]*mat[2][1]+mat[0][1]*mat[1][0]*mat[2][2]);


  double val1 = (mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2])/determinant;
  inverse_mat.push_back(val1);
  double val2 = -(mat[0][1]*mat[2][2]-mat[2][1]*mat[0][2])/determinant;
  inverse_mat.push_back(val2);
  double val3 = (mat[1][2]*mat[0][1]-mat[1][1]*mat[0][2])/determinant;
  inverse_mat.push_back(val3);
  double val4 = -(mat[2][2]*mat[1][0]-mat[2][0]*mat[1][2])/determinant;
  inverse_mat.push_back(val4);
  double val5 = (mat[0][0]*mat[2][2]-mat[2][0]*mat[0][2])/determinant;
  inverse_mat.push_back(val5);
  double val6 = -(mat[0][0]*mat[1][2]-mat[0][2]*mat[1][0])/determinant;
  inverse_mat.push_back(val6);
  double val7 = (mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0])/determinant;
  inverse_mat.push_back(val7);
  double val8 = -(mat[2][1]*mat[0][0]-mat[2][0]*mat[0][1])/determinant;
  inverse_mat.push_back(val8);
  double val9 = (mat[1][1]*mat[0][0]-mat[0][1]*mat[1][0])/determinant;
  inverse_mat.push_back(val9);


   cout << "INVERSE_MATRIX" << endl;
  for (int i = 0; i < 9; i++) cout <<inverse_mat[i] << " ";

  cout << " "<< endl;
  



  for (int i = 0; i < width; i++){
    for (int j = 0; j < height; j++) {
      double xA = best_H_matrix[0]*i + best_H_matrix[1]*j + best_H_matrix[2];
      double yA = best_H_matrix[3]*i + best_H_matrix[4]*j + best_H_matrix[5];
      double wA = best_H_matrix[6]*i + best_H_matrix[7]*j + best_H_matrix[8];


      double scaled_xA = xA/wA;
      double scaled_yA = yA/wA;
      // cout << scaled_xA << " "<< scaled_yA << endl;


      if (scaled_xA>0 && scaled_xA < width && scaled_yA >0 && scaled_yA<height){

          double red = (newA.Pixel(i,j).Red()+temp.Pixel(scaled_xA,scaled_yA).Red())/2;
          double green = (newA.Pixel(i,j).Green()+temp.Pixel(scaled_xA,scaled_yA).Green())/2;
          double blue = (newA.Pixel(i,j).Blue()+temp.Pixel(scaled_xA,scaled_yA).Blue())/2;

          R2Pixel *newPix = new R2Pixel(red,green,blue,1);
          SetPixel(i,j,*newPix); 


      } else {
          R2Pixel *newPix = new R2Pixel(1,1,1,1);
          SetPixel(i,j,*newPix); 

      }
    }
  }

  for (int a = 0 ; a < max_inlier_orig_set.size(); a++) {
    
       line(vecA[a][1],vecB[a][0],vecA[a][2],vecB[a][1],0,1,0);
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
    fprintf(stderr, "Unable to open image file: %s\n", filename);
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
    fprintf(stderr, "Unable to open image file: %s\n", filename);
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
    fprintf(stderr, "Unable to open image file: %s\n", filename);
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
      fprintf(stderr, "Unable to open image file: %s\n", filename);
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
      fprintf(stderr, "Unable to open image file: %s\n", filename);
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
    fprintf(stderr, "Unable to open image file: %s\n", filename);
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
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 100, TRUE);
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






