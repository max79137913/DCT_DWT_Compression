//*****************************************************************************
//
// Image.cpp : Defines the class operations on images
//
// Author - Parag Havaldar
// Code used by students as starter code to display and modify images
//
//*****************************************************************************

#include "Image.h"

#include "math.h"


#define pi 3.14159

int round(float x)
{
	int y = int (x+0.5);
	return y;
}

// Constructor and Desctructors
MyImage::MyImage() 
{
	Data = NULL;
	Width = -1;
	Height = -1;
	ImagePath[0] = 0;
}

MyImage::~MyImage()
{
	if ( Data )
		delete Data;
}


// Copy constructor
MyImage::MyImage( MyImage *otherImage)
{
	Height = otherImage->Height;
	Width  = otherImage->Width;
	Data   = new char[Width*Height*3];
	//strcpy(otherImage->ImagePath, ImagePath );
	strcpy( (char *)ImagePath, otherImage->ImagePath );
	for ( int i=0; i<(Height*Width*3); i++ )
	{
		Data[i]	= otherImage->Data[i];
	}


}



// = operator overload
MyImage & MyImage::operator= (const MyImage &otherImage)
{
	N=       otherImage.N;
	Height = otherImage.Height;
	Width  = otherImage.Width;
	Data   = new char[Width*Height*3];
	strcpy( (char *)ImagePath, otherImage.ImagePath );

	for ( int i=0; i<(Height*Width*3); i++ )
	{
		Data[i]	= otherImage.Data[i];
	}
	
	return *this;

}
	


// MyImage::ReadImage
// Function to read the image given a path
bool MyImage::ReadImage()
{

	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0 )
	{
		fprintf(stderr, "Image or Image properties not defined");
		fprintf(stderr, "Usage is `Image.exe Imagefile w h`");
		return false;
	}
	
	// Create a valid output file pointer
	FILE *IN_FILE;
	IN_FILE = fopen(ImagePath, "rb");
	if ( IN_FILE == NULL ) 
	{
		fprintf(stderr, "Error Opening File for Reading");
		return false;
	}

	// Create and populate RGB buffers
	int i;
	char *Rbuf = new char[Height*Width]; 
	char *Gbuf = new char[Height*Width]; 
	char *Bbuf = new char[Height*Width]; 

	for (i = 0; i < Width*Height; i ++)
	{
		Rbuf[i] = fgetc(IN_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		Gbuf[i] = fgetc(IN_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		Bbuf[i] = fgetc(IN_FILE);
	}
	
	// Allocate Data structure and copy
	Data = new char[Width*Height*3];
	for (i = 0; i < Height*Width; i++)
	{
		Data[3*i]	= Bbuf[i];
		Data[3*i+1]	= Gbuf[i];
		Data[3*i+2]	= Rbuf[i];
	}

	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(IN_FILE);

	return true;

}



// MyImage functions defined here
bool MyImage::WriteImage()
{
	// Verify ImagePath
	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0 )
	{
		fprintf(stderr, "Image or Image properties not defined");
		return false;
	}
	
	// Create a valid output file pointer
	FILE *OUT_FILE;
	OUT_FILE = fopen(ImagePath, "wb");
	if ( OUT_FILE == NULL ) 
	{
		fprintf(stderr, "Error Opening File for Writing");
		return false;
	}

	// Create and populate RGB buffers
	int i;
	char *Rbuf = new char[Height*Width]; 
	char *Gbuf = new char[Height*Width]; 
	char *Bbuf = new char[Height*Width]; 

	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = Data[3*i];
		Gbuf[i] = Data[3*i+1];
		Rbuf[i] = Data[3*i+2];
	}

	
	// Write data to file
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Rbuf[i], OUT_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Gbuf[i], OUT_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Bbuf[i], OUT_FILE);
	}
	
	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(OUT_FILE);

	return true;

}




// Here is where you would place your code to modify an image
// eg Filtering, Transformation, Cropping, etc.
bool MyImage::Modify()
{

	// TO DO by student
	
	// sample operation
	

	int i,j;
 // DCT transform
	for(j=0;j<Height;j=j+8)
	{
		for(i=0;i<Width;i=i+8)
		{
		
		//Dct(i,j);
		//DeDct(i,j);
		DctMatrix(i,j);
		}
	}

	//Dwt();
	//DctMatrix(0,0);
	//Dct(0,0);
	//DeDct(0,0);
	return false;
}
void MyImage::Dct(int startx , int starty )
{
	// Create and populate RGB buffers
	int i;
	int j;
	int *Rbuf = new int[Height*Width]; 
	int *Gbuf = new int[Height*Width]; 
	int *Bbuf = new int[Height*Width]; 



	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = int(Data[3*i]);
		Gbuf[i] = int(Data[3*i+1]);
		Rbuf[i] = int(Data[3*i+2]);
	}

	for (i = 0; i < Height*Width; i++)// convert to unsigned int
	{
		Bbuf[i] = Bbuf[i] &0xff;
		Gbuf[i] = Gbuf[i] &0xff; 
		Rbuf[i] = Rbuf[i] &0xff;
	}
	

	int u,v;

	float tmpr=0,tmpg=0,tmpb=0; 
	float Cu,Cv;

	for(v=starty;v<starty+8;v++){
			for(u=startx;u<startx+8;u++){
		
				
				if(u==0&&v==0){Cu=1/sqrt(2.0);Cv=1/sqrt(2.0);}
				else{Cu=1;Cv=1;}

				for(i=startx;i<startx+8;i++)
				{
					for(j=starty;j<starty+8;j++)
					{
					tmpr+=(float)Rbuf[i+j*Width]*cos(((2*i+1)*u*pi)/16)*cos(((2*j+1)*v*pi)/16);
					tmpg+=(float)Gbuf[i+j*Width]*cos(((2*i+1)*u*pi)/16)*cos(((2*j+1)*v*pi)/16);
					tmpb+=(float)Bbuf[i+j*Width]*cos(((2*i+1)*u*pi)/16)*cos(((2*j+1)*v*pi)/16);
					
					}
				}
								
				dctRbuf[u+v*Width]=(1.0/4.0)*Cu*Cv*tmpr;
				dctGbuf[u+v*Width]=(1.0/4.0)*Cu*Cv*tmpg;
				dctBbuf[u+v*Width]=(1.0/4.0)*Cu*Cv*tmpb;
	
				tmpr=tmpg=tmpb=0;

			}// u loop
	}//v loop



/*
	for (j = starty; j < starty+8; j++)// convert float to int
	{
		for(i=startx;i<startx+8;i++)
		{
		Bbuf[i+Width*j] = dctBbuf[i+Width*j];
		Gbuf[i+Width*j] = dctGbuf[i+Width*j]; 
		Rbuf[i+Width*j] = dctRbuf[i+Width*j];
		}
	}
	
	for (i = 0; i < Height*Width; i++)
	{
		Data[3*i]=char(Bbuf[i]);
		Data[3*i+1]=char(Gbuf[i]);
		Data[3*i+2]=char(Rbuf[i]);
	}
*/

	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
}
void MyImage::DeDct(int startx , int starty )
{
		int i,j;
		int u,v;

	float tmpr=0,tmpg=0,tmpb=0; 
	float Cu,Cv;

	int *Rbuf = new int[Height*Width]; 
	int *Gbuf = new int[Height*Width]; 
	int *Bbuf = new int[Height*Width]; 

	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = int(Data[3*i]);
		Gbuf[i] = int(Data[3*i+1]);
		Rbuf[i] = int(Data[3*i+2]);
	}

	

	for(j=starty;j<starty+8;j++){
			for(i=startx;i<startx+8;i++){


				
				

				for(u=startx;u<startx+8;u++)
				{
					for(v=starty;v<starty+8;v++)
					{
						
				if(u==0&&v==0){Cu=1/sqrt(2.0);Cv=1/sqrt(2.0);}
				else{Cu=1;Cv=1;}

				tmpr+=Cu*Cv*dctRbuf[u+v*Width]*cos(((2*i+1)*u*pi)/16)*cos(((2*j+1)*v*pi)/16);
				tmpg+=Cu*Cv*dctGbuf[u+v*Width]*cos(((2*i+1)*u*pi)/16)*cos(((2*j+1)*v*pi)/16);
				tmpb+=Cu*Cv*dctBbuf[u+v*Width]*cos(((2*i+1)*u*pi)/16)*cos(((2*j+1)*v*pi)/16);
					
					}
				}
				
				tmpr=(1.0/4.0)*tmpr;
				tmpg=(1.0/4.0)*tmpg;
				tmpb=(1.0/4.0)*tmpb;


				Bbuf[i+Width*j]=round(tmpb);
				Gbuf[i+Width*j]=round(tmpg);
				Rbuf[i+Width*j]=round(tmpr);
	
				tmpr=tmpg=tmpb=0;

				if(Bbuf[i+Width*j]<0){Bbuf[i+Width*j]=0;}
				else if(Bbuf[i+Width*j]>255){Bbuf[i+Width*j]=255;}
				if(Gbuf[i+Width*j]<0){Gbuf[i+Width*j]=0;}
				else if(Gbuf[i+Width*j]>255){Gbuf[i+Width*j]=255;}
				if(Rbuf[i+Width*j]<0){Rbuf[i+Width*j]=0;}
				else if(Rbuf[i+Width*j]>255){Rbuf[i+Width*j]=255;}

				

			}// i loop
	}//j loop


	
	for (i = 0; i < Height*Width; i++)
	{
		




		Data[3*i]=char(Bbuf[i]);
		Data[3*i+1]=char(Gbuf[i]);
		Data[3*i+2]=char(Rbuf[i]);
	}

	delete Rbuf;
	delete Gbuf;
	delete Bbuf;

}
void MyImage::DctMatrix(int startx, int starty) // Dct encode
{
	float U[8][8];
	float UT[8][8];
	float A[3][8][8];
	int x,y;

////////////////////////////// init RGB
	int *Rbuf = new int[Height*Width]; 
	int *Gbuf = new int[Height*Width]; 
	int *Bbuf = new int[Height*Width]; 

	int i;
	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = int(Data[3*i]);
		Gbuf[i] = int(Data[3*i+1]);
		Rbuf[i] = int(Data[3*i+2]);
	}
		for (i = 0; i < Height*Width; i++)// convert to unsigned int
	{
		Bbuf[i] = Bbuf[i] &0xff;
		Gbuf[i] = Gbuf[i] &0xff; 
		Rbuf[i] = Rbuf[i] &0xff;
	}
//////////////////////////////////// making U matrix
	for(y=0;y<8;y++)
	{
		for(x=0;x<8;x++)
		{
			float value;
			if(y==0){value=sqrt(2.0)/4;}
			else
			{
			value=cos(((y+2*x*y)*pi)/16)/2;		
			}

		U[x][y]=value;
			
		}
	}
/////////////////////////////////making UT matrix
 
	for(y=0;y<8;y++)
	{
		for(x=0;x<8;x++)
		{
			float value;
			if(x==0){value=sqrt(2.0)/4;}
			else
			{
			value=cos(((x+2*x*y)*pi)/16)/2;		
			}

		UT[x][y]=value;
			
		}
	}
/////////////////////////////////
	i=0;
	int j=0;
	for(y=starty;y<starty+8;y++)
	{
		
		for(x=startx;x<startx+8;x++)
		{

	A[0][i][j]=(float)Bbuf[x+y*Width];
	A[1][i][j]=(float)Gbuf[x+y*Width];
	A[2][i][j]=(float)Rbuf[x+y*Width];
		i++;
		}
		i=0;
		j++;
	}

	Matrix EUendmatrix[3];// encode
	Matrix EUTendmatrix[3];

	Matrix DUendmatrix[3];//decode
	Matrix DUTendmatrix[3];
	/////////////////////// A[0]->B
	ComputeMatrix(U,A[0],&EUendmatrix[0]);
	ComputeMatrix(EUendmatrix[0],UT,&EUTendmatrix[0]);
	///////////////////////////
	DctZig(&EUTendmatrix[0]);  //B->zigB
	////////////////////////////zigB->A[0]
	ComputeMatrix(UT,EUTendmatrix[0],&DUTendmatrix[0]);
	ComputeMatrix(DUTendmatrix[0],U,&DUendmatrix[0]);
	//////////////////////////
	/////////////////////// A[1]->B
	ComputeMatrix(U,A[1],&EUendmatrix[1]);
	ComputeMatrix(EUendmatrix[1],UT,&EUTendmatrix[1]);
	///////////////////////////
	DctZig(&EUTendmatrix[1]);  //B->zigB
	////////////////////////////zigB->A[1]
	ComputeMatrix(UT,EUTendmatrix[1],&DUTendmatrix[1]);
	ComputeMatrix(DUTendmatrix[1],U,&DUendmatrix[1]);
	////////////////////////////
	/////////////////////// A[2]->B
	ComputeMatrix(U,A[2],&EUendmatrix[2]);
	ComputeMatrix(EUendmatrix[2],UT,&EUTendmatrix[2]);
	///////////////////////////
	DctZig(&EUTendmatrix[2]);  //B->zigB
	////////////////////////////zigB->A[2]
	ComputeMatrix(UT,EUTendmatrix[2],&DUTendmatrix[2]);
	ComputeMatrix(DUTendmatrix[2],U,&DUendmatrix[2]);
	////////////////////////////
	
	i=0;j=0;// back to 0
		for(y=starty;y<starty+8;y++)
	{
		
		for(x=startx;x<startx+8;x++)
		{

	Bbuf[x+y*Width]=round(DUendmatrix[0][i][j]);
	Gbuf[x+y*Width]=round(DUendmatrix[1][i][j]);
	Rbuf[x+y*Width]=round(DUendmatrix[2][i][j]);

				if(Bbuf[x+y*Width]<0){Bbuf[x+y*Width]=0;}
				else if(Bbuf[x+y*Width]>255){Bbuf[x+y*Width]=255;}
				if(Gbuf[x+y*Width]<0){Gbuf[x+y*Width]=0;}
				else if(Gbuf[x+y*Width]>255){Gbuf[x+y*Width]=255;}
				if(Rbuf[x+y*Width]<0){Rbuf[x+y*Width]=0;}
				else if(Rbuf[x+y*Width]>255){Rbuf[x+y*Width]=255;}
		i++;
		}
		i=0;
		j++;
	}
		i=0;j=0;
		for (i = 0; i < Height*Width; i++)
	{
		Data[3*i]=char(Bbuf[i]);
		Data[3*i+1]=char(Gbuf[i]);
		Data[3*i+2]=char(Rbuf[i]);
	}

	delete Rbuf;
	delete Gbuf;
	delete Bbuf;


}
void MyImage::ComputeMatrix(Matrix first, Matrix second, Matrix *EndMatrix)
{
	
	
			Matrix FMatrix;
			int x,y;
	
		for(y=0;y<8;y++){  //multify push matrix

		for(x=0;x<8;x++)
	{
		FMatrix[x][y]=first[0][y]*second[x][0]+first[1][y]*second[x][1]+first[2][y]*second[x][2]+first[3][y]*second[x][3]+first[4][y]*second[x][4]+first[5][y]*second[x][5]+first[6][y]*second[x][6]+first[7][y]*second[x][7];
	}		
		
					    }
/////////////////////////////////////////////////
	for(y=0;y<8;y++){   /// assign F to EndMatrix

			for(x=0;x<8;x++)
			{
			(*EndMatrix)[x][y]=FMatrix[x][y];
			}
											}

}
void MyImage::DctZig(Matrix *zig)
{

	// fill along diagonal stripes (oriented as "/")
	int dimension=8;
	int lastValue = 8 * 8 - 1;
	int currDiag = 0;
	int From;
	int To;
	int i;
	int row;
	int col;
	int times=0;

	int m;
	m=round(N/4096);


	do
	{
		if ( currDiag < dimension ) // if doing the upper-left triangular half
		{
			From = 0;
			To = currDiag;
		}
		else // doing the bottom-right triangular half
		{
			From = currDiag - dimension + 1;
			To = dimension - 1;
		}
 
		for ( i = From; i <= To; i++ )
		{
			if ( currDiag % 2 == 0 ) // want to fill upwards
			{
				col = To - i + From;
				row = i;
			}
			else // want to fill downwards
			{
				col = i;
				row = To - i + From;
			}
 
			if(times>m-1)
			{(*zig)[row][col] =0;}
			times++;
		}
		
		currDiag++;
	}
	while ( currDiag <= 14 );
}
void MyImage::Dwt()
{

	int counter;
	////////////////////////////// init RGB
	int *Rbuf = new int[Height*Width]; 
	int *Gbuf = new int[Height*Width]; 
	int *Bbuf = new int[Height*Width]; 

	float *Btmp = new float[Height*Width]; 
	float *Gtmp = new float[Height*Width]; 
	float *Rtmp = new float[Height*Width]; 

	float *Binter = new float[Height*Width]; 
	float *Ginter = new float[Height*Width]; 
	float *Rinter = new float[Height*Width]; 
    

	int i,j,k;
	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = int(Data[3*i]);
		Gbuf[i] = int(Data[3*i+1]);
		Rbuf[i] = int(Data[3*i+2]);
	}
		for (i = 0; i < Height*Width; i++)// convert to unsigned int
	{
		Bbuf[i] = Bbuf[i] &0xff;
		Gbuf[i] = Gbuf[i] &0xff; 
		Rbuf[i] = Rbuf[i] &0xff;
	}

///////////////////////////////////////

		for (i = 0; i < Height*Width; i++)
	{
		Binter[i]=(float)Bbuf[i];
		Ginter[i]=(float)Gbuf[i];
		Rinter[i]=(float)Rbuf[i];
	}

	int x=0,y=0; 
	int scale=Width;



/////////////////////////////////////////row transform
	for(k=0;k<9;k++)
	{
		for(j=0;j<Height;j++)
		{
			for(i=0;i<scale;i=i+2)
			{		

			//Btmp[Width/2+x][j]=(float)Bbuf[i+j*Width]-(float(Bbuf[i+j*Width])+float(Bbuf[(i+1)+j*Width]))/2;
			Btmp[x+y*Width]=(Binter[i+j*Width]+Binter[(i+1)+j*Width])/2;
			Btmp[(scale/2+x)+y*Width]=Binter[i+j*Width]-Btmp[x+y*Width];
			//////////////////////////////////////////////////////////////////////
			Gtmp[x+y*Width]=(Ginter[i+j*Width]+Ginter[(i+1)+j*Width])/2;
			Gtmp[(scale/2+x)+y*Width]=Ginter[i+j*Width]-Gtmp[x+y*Width];
			//////////////////////////////////////////////////////////////////////
			Rtmp[x+y*Width]=(Rinter[i+j*Width]+Rinter[(i+1)+j*Width])/2;
			Rtmp[(scale/2+x)+y*Width]=Rinter[i+j*Width]-Rtmp[x+y*Width];
			x++;
			}	
			x=0;
			y++;
		}
		y=0;
		scale=scale/2;
		////////////////FINISH row transform onetime

		for (counter = 0; counter < Height*Width; counter++)
	{
		Binter[counter]=Btmp[counter];
		Ginter[counter]=Gtmp[counter];
		Rinter[counter]=Rtmp[counter];
		
	}


	}
///////////////////////////////////////
		x=0;y=0; //end row dwt
///////////////////////////////////////

for (counter = 0; counter < Height*Width; counter++)
	{
		Binter[counter]=Btmp[counter];
		Ginter[counter]=Gtmp[counter];
		Rinter[counter]=Rtmp[counter];
		
	}

/////////////////////////////////////////Dwt column transform
scale=Height;
	

	for(k=0;k<9;k++)
	{
		for(i=0;i<Width;i++)
		{
			for(j=0;j<scale;j=j+2)
			{	

			//Btmp[Width/2+x][j]=(float)Bbuf[i+j*Width]-(float(Bbuf[i+j*Width])+float(Bbuf[(i+1)+j*Width]))/2;
			Btmp[x+y*Width]=(Binter[i+j*Width]+Binter[(i)+(j+1)*Width])/2;
			Btmp[(x)+(y+scale/2)*Width]=Binter[i+j*Width]-Btmp[x+y*Width];
			//////////////////////////////////////////////////////////////////////
			Gtmp[x+y*Width]=(Ginter[i+j*Width]+Ginter[(i)+(j+1)*Width])/2;
			Gtmp[(x)+(y+scale/2)*Width]=Ginter[i+j*Width]-Gtmp[x+y*Width];
			//////////////////////////////////////////////////////////////////////
			Rtmp[x+y*Width]=(Rinter[i+j*Width]+Rinter[(i)+(j+1)*Width])/2;
			Rtmp[(x)+(y+scale/2)*Width]=Rinter[i+j*Width]-Rtmp[x+y*Width];

			y++;
			}	
			y=0;
			x++;
		}
		x=0;
		scale=scale/2;
		////////////////FINISH column transform onetime

		for (counter = 0; counter < Height*Width; counter++)
	{
		Binter[counter]=Btmp[counter];
		Ginter[counter]=Gtmp[counter];
		Rinter[counter]=Rtmp[counter];
		
	}


	}
////////////////////////////////////
		x=0;y=0; //end column dwt
////////////////////////////////////

////////////////////////////////////// dwt zigzag



	// fill along diagonal stripes (oriented as "/")
	int dimension=512;
	int lastValue = 512 * 512 - 1;
	int currDiag = 0;
	int From;
	int To;
	int row;
	int col;
	int times=0;

	int m;
	m=N;


	do
	{
		if ( currDiag < dimension ) // if doing the upper-left triangular half
		{
			From = 0;
			To = currDiag;
		}
		else // doing the bottom-right triangular half
		{
			From = currDiag - dimension + 1;
			To = dimension - 1;
		}
 
		for ( i = From; i <= To; i++ )
		{
			if ( currDiag % 2 == 0 ) // want to fill upwards
			{
				col = To - i + From;
				row = i;
			}
			else // want to fill downwards
			{
				col = i;
				row = To - i + From;
			}
 
			if(times>m-1)
			{
				Btmp[row+col*Width]=0;
				Gtmp[row+col*Width]=0;
				Rtmp[row+col*Width]=0;
			}
			times++;
		}
		
		currDiag++;
	}
	while ( currDiag <= 1022 );



//////////////////////////////////////

//////////////////////////////////////    IDwt Column
		for (i = 0; i < Height*Width; i++)
	{
		Binter[i]=Btmp[i];
		Ginter[i]=Gtmp[i];
		Rinter[i]=Rtmp[i];
		
	}


/////////////////////////////////////// 

scale=1;

for(i=0;i<Width;i++)
{
	for(k=0;k<9;k++)//level
	{
		for(j=0;j<scale;j++)
		{
		Btmp[(i)+(2*j)*Width]=Binter[i+j*Width]+Binter[(i)+(j+scale)*Width];
		Btmp[(i)+(2*j+1)*Width]=Binter[i+j*Width]-Binter[(i)+(j+scale)*Width];
		///////////////////////////////////////////////////////////////
		Gtmp[(i)+(2*j)*Width]=Ginter[i+j*Width]+Ginter[(i)+(j+scale)*Width];
		Gtmp[(i)+(2*j+1)*Width]=Ginter[i+j*Width]-Ginter[(i)+(j+scale)*Width];
		///////////////////////////////////////////////////////////////
		Rtmp[(i)+(2*j)*Width]=Rinter[i+j*Width]+Rinter[(i)+(j+scale)*Width];
		Rtmp[(i)+(2*j+1)*Width]=Rinter[i+j*Width]-Rinter[(i)+(j+scale)*Width];
		
		}///j loop

		scale=2*scale;
		for (counter = 0; counter < Height*Width; counter++)// save content
		{
		Binter[counter]=Btmp[counter];
		Ginter[counter]=Gtmp[counter];
		Rinter[counter]=Rtmp[counter];
		
		}
	}
	scale=1;

}
////////////////////////////// end IDWT Column



//////////////////////////////////////    IDwt Row
		for (i = 0; i < Height*Width; i++)
	{
		Binter[i]=Btmp[i];
		Ginter[i]=Gtmp[i];
		Rinter[i]=Rtmp[i];
		
	}

/////////////////////////////////////// 

scale=1;

for(j=0;j<Height;j++)
{
	for(k=0;k<9;k++)//level
	{
		for(i=0;i<scale;i++)
		{
		Btmp[(2*i)+j*Width]=Binter[i+j*Width]+Binter[(i+scale)+j*Width];
		Btmp[(2*i+1)+j*Width]=Binter[i+j*Width]-Binter[(i+scale)+j*Width];
		///////////////////////////////////////////////////////////////
		Gtmp[(2*i)+j*Width]=Ginter[i+j*Width]+Ginter[(i+scale)+j*Width];
		Gtmp[(2*i+1)+j*Width]=Ginter[i+j*Width]-Ginter[(i+scale)+j*Width];
		///////////////////////////////////////////////////////////////
		Rtmp[(2*i)+j*Width]=Rinter[i+j*Width]+Rinter[(i+scale)+j*Width];
		Rtmp[(2*i+1)+j*Width]=Rinter[i+j*Width]-Rinter[(i+scale)+j*Width];
		
		}///i loop

		scale=2*scale;
		for (counter = 0; counter < Height*Width; counter++)// save content
		{
		Binter[counter]=Btmp[counter];
		Ginter[counter]=Gtmp[counter];
		Rinter[counter]=Rtmp[counter];
		
		}
	}
	scale=1;

}
////////////////////////////// end IDWT Row


	/////////////store into Data
	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i]=round(Btmp[i]);
		Gbuf[i]=round(Gtmp[i]);
		Rbuf[i]=round(Rtmp[i]);



		if(Bbuf[i]<0){Bbuf[i]=0;}
		else if(Bbuf[i]>255){Bbuf[i]=255;}
		if(Gbuf[i]<0){Gbuf[i]=0;}
		else if(Gbuf[i]>255){Gbuf[i]=255;}
		if(Rbuf[i]<0){Rbuf[i]=0;}
		else if(Rbuf[i]>255){Rbuf[i]=255;}
		
	}

/*
		for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i]=int(Btmp[i]);
		Gbuf[i]=int(Gtmp[i]);
		Rbuf[i]=int(Rtmp[i]);
		
	}
*/
	for (i = 0; i < Height*Width; i++)
	{
		Data[3*i]=char(Bbuf[i]);
		Data[3*i+1]=char(Gbuf[i]);
		Data[3*i+2]=char(Rbuf[i]);
	}

	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	delete Rtmp;
	delete Gtmp;
	delete Btmp;
	delete Rinter;
	delete Ginter;
	delete Binter;

}
