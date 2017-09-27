// QRTestSuite.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <chrono>
#include <zbar.h>
using namespace std;
using namespace cv;
using namespace zbar;
//###############################################################################
// Initialization
//###############################################################################

//#################		Tags
const bool benchmark = true;
const String benchmarktype = "stream"; //stream mp4, jpg, png <--- prob useless since no gain
const bool debug = true;
const bool showCalc = true;

//#################		Global Variables
Mat image; //back to multiFIPdetector
int64 e1, e2;
float t;
vector<float> allFPS;
float avgFPS = 0.0, highFPS = 0.0, lowFPS = 999.0, globalFPS = 0.0; //FPS Performance Calculation happening global so minimal impact on algorithm
vector<int> fips; //just saves the number of FiPs detected for Benchmark
Point prevPos;

//#################		Classes
class FiP;
class QRCode;

//#################		Main Methods
int QRTest();
int fpsCalc(float fps);
//int imageLoadTest();

//Detector Methods
vector<FiP> cv_FiPdetection(Mat inputImage, vector<FiP> prevImage);
QRCode cv_QRdetection(vector<FiP> fipImage, QRCode qrpPrevImage);
int cv_HarrisCorner(Mat img);
vector<FiP> cv_getFiPOrder(vector<FiP> unordered);

//#################		Support Methods
int cv_findCorners(Point& pA, FiP fip_B, FiP fip_C, Point& pD, Point QRPos);
Point cv_getOuterCorner(FiP fip, Point center);
Point cv_getOuterCorner(FiP fip, Point center, int& index);
bool cv_getIntersection(Point a1, Point a2, Point b1, Point b2, Point& r);
float cv_lineLineAngle(Point l1_1, Point l1_2, Point l2_1, Point l2_2);
float cv_vectorSize(Point a);
Point cv_getCentroid(vector<Point> contour);
float cv_euclideanDist(Point p, Point q);
bool cv_inRect(vector<Point> rectangle, Point p);
bool cv_inFiPRegTesting(vector<vector<Point>> &FiPRegs, vector<Point> Contour, vector<bool> &updated);
int cv_CandidateInRegion(vector<Point> contour, vector<vector<Point> > candidates);
int cv_outputHisto(Mat input);
float cross(Point2f v1, Point2f v2);
String cv_getOrientation(Point a, Point b);
//Decode Function
int decode(Mat inputImage);

//###############################################################################
//Methods
//###############################################################################
class FiP {
public:
	FiP();
	FiP (int pRelPos);
	FiP (Point pPos, vector<Point> pShape);
	Point pos;
	vector<Point> shape;
	int relPos; //relative pos in the QR code: 0 the corner, 1 and 2 the diagonal Fips: unkown -1? or 0?

	void setRelPos(int pRelPos) {
		relPos = pRelPos;
	}

};

FiP::FiP() {
	
}

FiP::FiP(int pRelPos) {
	relPos = pRelPos;
}

FiP::FiP (Point pPos, vector<Point> pShape) {
	pos = pPos;
	shape = pShape;
}

class QRCode { //Class for the QR-Code that saves all relevant information for one Instance
	public:
		QRCode(int ptag);
		Point pos; // x and y pos of the middle of the QR-Code (mostly for testing as the Corners will be used for processing)
		int tag; // Tag that will identify this QR-Code so that every QR-data will only have one corresponding QR-Code
		FiP fip_a, fip_b, fip_c; //The finderpatterns
		Point a, b, c, d; //the Corner Points
		String orientation; //does it make sense to create a new class for this?
		//int screenSize; //Size in ScreenSpace
		//float realSize; //Size in Centimeters - probably 20cm 

};

QRCode::QRCode(int ptag) {
	tag = ptag;
}

int main() {
	for (int i = 0; i < 10; i++)
	{	
		int out = QRTest();
		if (out == -1) {
			cerr << "ERR: Manual Abort" << endl;
			break;
		}
		fpsCalc(globalFPS);
	}
	cout << "--------------Complete Run finished--------------" << endl;
	cout << "Average FPS: " << avgFPS << endl;
	cout << "Highes FPS: " << highFPS << endl;
	cout << "Lowest FPS: " << lowFPS << endl;
	waitKey(0);
}

int fpsCalc(float fps) {
	float sumFPS = 0.0;
	allFPS.push_back(fps);

	for (int i = 0; i < allFPS.size(); i++)
		sumFPS += allFPS.at(i);
	avgFPS = sumFPS / allFPS.size();

	if (fps > highFPS)
		highFPS = fps;
	else if (fps < lowFPS)
		lowFPS = fps;

	return 0;
}

int QRTest()
{
	e1 = getTickCount();
	cout << "--------------Test Start--------------" << endl;
	//Variable Initializations
	int key = 0; //input key of the hughgui classbb
	VideoCapture capture;//("C:/Users/Frederik/Desktop/VidTests/720-frame-jpg/image-%05d.jpg");
	vector<int> timePassedAll;
	vector<FiP> FiPList;
	QRCode QRList = QRCode(0);
	int solo = 0, partial = 0, full = 0, over = 0;
	int frame = 0;

	// Creation of Intermediate 'Image' Objects required later
	Mat empty(Size(100, 100), CV_MAKETYPE(image.depth(), 1));

	if (benchmark == false) {
		capture.open(0);
	}
	else {
		if (benchmarktype == "stream")
			capture.open("C:/Users/Frederik/Desktop/VidTests/long-exposure.mp4");
			//capture.open("C:/Users/Frederik/Desktop/VidTests/720p-dist-move.mp4");
		else if (benchmarktype == "jpg")
			capture.open("C:/Users/Frederik/Desktop/VidTests/720-frame-jpg/image-%05d.jpg");
		else if (benchmarktype == "png")
			capture.open("C:/Users/Frederik/Desktop/VidTests/720-frame-png/image-%07d.png");
	}

	if (!capture.isOpened()) {
		cerr << " ERR: Unable find input Video source." << endl;
		return -1;
	}

	while (key != 'q' && capture.isOpened()) //Main Computation Loop
	{
		frame++; 

		//Pause for Observing Output - This will disrupt the timing
		if (key == 'p') {
			key = waitKey(0);
		}

		capture >> image;

		if (image.empty()) { // If there is no more pictures in the benchmark the loop will break;
			//Benchmark  Time
			e2 = getTickCount();
			t = (e2 - e1) / getTickFrequency();
			int frames = capture.get(CAP_PROP_FRAME_COUNT);
			float fps = (frames / t);
			globalFPS = fps;
			//Benchmark Output Performance

			//Output
			cout << "-------End of Videofile reached-------" << endl;
			cout << "Frames Computed       : " << frames << endl;
			cout << "-------------Time Metrics-------------" << endl;
			cout << "Time Passed           : " << t << endl;
			cout << "Frames per Sec        : " << fps << endl;
			cout << "-----------Accuracy Metrics-----------" << endl;
			cout << "Frames solo           : " << ((float)solo / (float)frames) * 100 << "%" << endl;  //one found
			cout << "Frames partial        : " << ((float)partial / (float)frames) * 100 << "%" <<endl;  //two found
			cout << "Frames full detection : " << ((float)full / (float)frames) * 100 << "%" << endl; //Three FiPs corrected or Reconstructed
			cout << "Frames overshoot      : " << over << endl;
			//precicion?? <- we need ground truth
			//QR Codes?? <- For multiple
			//somehow measure how correct the guesses were ~ similar to precicion? 

			return 0;
		}

		//Reset of Variables
		//vector<vector<Point> > finderPatterns;
		//vector<vector<Point> > finderCandidate;


		//Image Processing
		FiPList = cv_FiPdetection(image, FiPList);
		QRList = cv_QRdetection(FiPList, QRList);



		
		int fipAmount = FiPList.size();
		fips.push_back(fipAmount);
		if (fipAmount == 3) {
			full++;
			partial++;
			solo++;
		}
		else if (fipAmount == 2) {
			partial++;
			solo++;
		}
		else if (fipAmount == 1) {
			solo++;
		}
		if (fipAmount > 3) {
			over++;
		}

		//Debug Rendering
		if (showCalc)
			imshow("image", image);
		else
			imshow("image", empty);

		key = waitKey(1);	// OPENCV: wait for 1ms before accessing next frame
	} //End of Main Computation Loop
	if (key == 'q')
		return -1;
	else
		return 0;
}

/* calling a method an reinitializing a lot of Variables should create a overhead compared to just use one method?
	Maybe I should do a speed test since it would be far easier to understand in more methods but if it costs more time
	nothing is really gained*/
vector<FiP> cv_FiPdetection(Mat inputImage, vector<FiP> prevImage) /*
									input:

									output:

									*/
{
	vector<FiP> fipList;
	vector<vector<Point> > contours;
	vector<vector<Point> > fiPReg;
	vector<Vec4i> hierarchy;
	vector<Point> pointsseq;    //used to save the approximated sides of each contour
	vector<vector<Point> > finderPatterns;
	vector<vector<Point> > finderCandidate;
	vector<bool> updated;
	updated.clear();

	cvtColor(inputImage, inputImage, CV_RGB2GRAY);										// Convert Image captured from Image Input to GrayScale	
	
	resize(inputImage, inputImage, Size(inputImage.size().width / 2, inputImage.size().height / 2));
	threshold(inputImage, inputImage, 100, 150, THRESH_BINARY);							//<- Probably some kind of local threshold better in the Areas of Interest 
	resize(inputImage, inputImage, Size(inputImage.size().width * 2, inputImage.size().height * 2));

	//imshow("debug", inputImage);

	Canny(inputImage, inputImage, 80, 150, 3);											// Apply Canny edge detection on the gray image
	findContours(inputImage, contours, hierarchy, RETR_TREE, CHAIN_APPROX_TC89_KCOS);



	for (int i = 0; i < contours.size(); i++)
	{

		//Find the approximated polygon of the contour we are examining
		approxPolyDP(contours[i], pointsseq, arcLength(contours[i], true)*0.02, true); // <- this somehow doesnt work
		if (pointsseq.size() == 4)      // only quadrilaterals contours are examined // <- this neither or I misunderstand it
		{
			int k = i;
			int c = 0;

			while (hierarchy[k][2] != -1)
			{
				//cout << "RUNS :  " << testo << "\n";
				k = hierarchy[k][2];
				c = c + 1;
				//testo++;
			}
			/*if (hierarchy[k][2] != -1)
			cout << "WARNING";
			c = c + 1;*/

			if (c >= 5)
			{
				finderPatterns.push_back(contours[k]);
				/*
				if (mark == 0)		A = i;
				else if (mark == 1)	B = i;		// i.e., A is already found, assign current contour to B
				else if (mark == 2)	C = i;		// i.e., A and B are already found, assign current contour to C
				mark = mark + 1;*/
				//cout << "FiP area:  " << contourArea(contours[k]) << "\n";

				//if Center not in ANY FiPreg Add to FiPReg
				//if (!cv_inFiPReg(fiPReg, cv_getCentroid(contours[k - c]))) { // <-- will try to do the update in here? or change it do give back the index of the fitting rectangle for later use
				// What if more Points are included in a region or a Point in more then one? second shouldn't really happen 
				// the first part is only for the Candidates
				vector<Point> fipSquare;
				approxPolyDP(contours[k - c], fipSquare, arcLength(contours[i], true)*0.02, true);
				if (!cv_inFiPRegTesting(fiPReg, fipSquare, updated)) {
					fiPReg.push_back(fipSquare);
					//updated[updated.size() - 1].flip();// = true;
					updated.push_back(true);
					//cout << fiPReg.size() << "\n";
				}



				//if Center in a FiPReg //TODO: decide if marked as found or put a connection so that possible multiples can be detected but that shouldnt happen so how would the problem be solved?
			}
			if (c = 4) {
				approxPolyDP(contours[k], pointsseq, arcLength(contours[i], true)*0.02, true);
				if (pointsseq.size() == 4) {
					if (isContourConvex(pointsseq)) {
						if (contourArea(contours[k]) > 10) { //TODO: check how expensive and finetune value? find less magic number
							finderCandidate.push_back(contours[k]);
							//cout << "Candidate area:  " << contourArea(contours[k]) << "\n";
							if (showCalc == true) {
								//drawContours(image, vector<vector<Point> >(1, pointsseq), -1, Scalar(0, 255, 0), 1, 8);
							}
							//size++;
						}
						//convex++;
					}
					//points++;
				}
				//all++;


			}

		}
	}

	for (auto &fip : fiPReg) {
		Point Center = cv_getCentroid(fip);
		fipList.push_back(FiP(Center, fip));
	}

	return fipList;
}

QRCode cv_QRdetection(vector<FiP> fipImage, QRCode qrPrevImage) {
	Point QRPos;
	FiP fip_A, fip_B, fip_C;
	bool found = false;
	// process FiPs
	fipImage = cv_getFiPOrder(fipImage);
	// if 3 or 2 diagonal we already have the position
	if (fipImage.size() == 3) {
		if (fipImage[0].relPos == 0) {
			QRPos = (fipImage[1].pos + fipImage[2].pos)*.5;
			fip_A = fipImage[0];
			fip_B = fipImage[1];
			fip_C = fipImage[2];
		}
		else if (fipImage[1].relPos == 0) {
			QRPos = (fipImage[0].pos + fipImage[2].pos)*.5;
			fip_A = fipImage[1];
			fip_B = fipImage[0];
			fip_C = fipImage[2];
		}
		else if (fipImage[2].relPos == 0) {
			QRPos = (fipImage[0].pos + fipImage[1].pos)*.5;
			fip_A = fipImage[1];
			fip_B = fipImage[2];
			fip_C = fipImage[0];
		}
		found = true;
	}
	else if (fipImage.size() == 2) {
		if (fipImage[0].relPos != 0 && fipImage[1].relPos != 0) {
			QRPos = (fipImage[0].pos + fipImage[1].pos)*.5;
			fip_A = FiP(-1);
			fip_B = fipImage[0];
			fip_C = fipImage[1];
			found = true;
		}
	}
	// if 2 parallel possible position
	// What could we do?
	// check if there is a tag expected in one of the two possible position and then accept that
	// if 1 check possible positions
	// What could we do?
	// check if there is a tag expected in the radius?
	// ignore might work best? (only if this happens rarely which seems to be the case)
	// check if there is tag for the Pos



	// if yes load data from app

	// if no reconstruct planar QR code and load data

	// check if tag exists else give new one
	// we want a distance to a tag from a previous image to check if this is the same tag
	//	this could create problems when there are overlaying tags
	// It would be clever to have this distance to change depending on the movement speed/ turning speed
	// With a normal movement the distance between the centers rarely goes over 50 pixels but with rapid movements it tends to get far higher
	// can we approximate the speed just by the picture in a really fast way? Or should we just describe that here a injection of attributes send by
	// the phone would be extremly helpfull
	// Because actually the detection of same Tag will be MOST Important when fast movements take place -> distortion of the picture
	// we could see how other detected FiPs or QRs behave but we still would have to have method to make sure which FiP moved where espacially in the
	// case that some might disappear 
	// Also in the case of multi code we could tag all untagged QRs simply so that all did a similar movement -> we say that a quick forward/backward
	// movement is not as likely as a rotational movement

	if (!found) { //This is for now since there is no reconstruction so only if I can actually do the calculations we continue
		//cout << "not Found" << endl;
		return qrPrevImage; 
	}
	//if (cv_vectorSize(qrPrevImage.pos - QRPos) <= 150.0) {
	if ((cv_euclideanDist(qrPrevImage.pos, QRPos)) <= 150.0) {
		//cout << "direct" << endl;
		return qrPrevImage;
	}
	else {
		//Identify
		
		//get CornerPointA in Case it exists already
		Point pA;
		if (fip_A.relPos == -1)
			pA = Point(0,0);
		else 
			pA = cv_getOuterCorner(fip_A, QRPos);
		//get Corners A if it doesnt and D in any Case
		Point pD;
		int cornersOut = cv_findCorners(pA, fip_B, fip_C, pD, QRPos);
		//cout << "Calculated Corners" << endl;
		//Now with the corners Reconstruct the planar image of the QRCode

		//TODO This only works with Point2f which isnt Paintable also our QR Image is quite Poor maybe we should get it from the original image 
		//and not from the changed one to find the FIP -- but here we should see first what zBar or another library can make out of it 
		vector<Point2f> qrImg, dst;
		qrImg.push_back(Point2f(pA));
		qrImg.push_back(Point2f(cv_getOuterCorner(fip_B,  QRPos)));
		qrImg.push_back(Point2f(pD)); //<- order important?
		qrImg.push_back(Point2f(cv_getOuterCorner(fip_C, QRPos)));

		//For Testing Do different later!
		Mat qr, qr_raw, qr_gray, qr_gray_sharp, qr_thres;
		qr_raw = Mat::zeros(200, 200, CV_8UC3);
		qr = Mat::zeros(200, 200, CV_8UC3);
		qr_gray = Mat::zeros(200, 200, CV_8UC1);
		qr_gray_sharp = Mat::zeros(200, 200, CV_8UC1);
		qr_thres = Mat::zeros(200, 200, CV_8UC1);

		//replace with non magic Numbers?
		dst.push_back(Point2f(0, 0));
		dst.push_back(Point2f(qr.cols, 0));
		dst.push_back(Point2f(qr.cols, qr.rows));
		dst.push_back(Point2f(0, qr.rows));

		Mat warp_matrix;

		//cout << "draw" << endl;
		//cout << qrImg.size() << endl;
		//cout << dst.size() << endl;
		//drawContours(image, vector<vector<Point2f> >(1.0, qrImg), -1, Scalar(0, 0, 255), 1, 8);
		//imshow("contour", image);
		//waitKey(0);
		if (qrImg.size() == 4 && dst.size() == 4)			// Failsafe for WarpMatrix Calculation to have only 4 Points with src and dst
		{
			warp_matrix = getPerspectiveTransform(qrImg, dst);
			warpPerspective(image, qr_raw, warp_matrix, Size(200, 200));
			copyMakeBorder(qr_raw, qr, 10, 10, 10, 10, BORDER_CONSTANT, Scalar(255, 255, 255));

			cvtColor(qr, qr_gray, CV_RGB2GRAY);
			//sharpen Image attempt Gaussian Blurr failed
			//GaussianBlur(qr_gray, qr_gray_sharp, Size(3,3), 1);
			//addWeighted(qr_gray, 1.5, qr_gray_sharp, -0.5, 0, qr_gray_sharp);
			//threshold(qr_gray_sharp, qr_thres, 130, 255, CV_THRESH_BINARY);
			//imshow("QR code", qr_thres);
			//check image quality before processing?
			threshold(qr_gray, qr_thres, 127, 255, CV_THRESH_BINARY);
			imshow("QR code old", qr_thres);
			//waitKey(0);
			decode(qr_thres);
		}




		//If tag exists retag
		int tag = qrPrevImage.tag + 1; //<- usually find next free tag - global variable?
		QRCode newCode(tag);
		newCode.pos = QRPos;
		if (fipImage.size() == 3)
			newCode.orientation = cv_getOrientation(pA, pD); //<- only do this if we had all threeFiPs so its sure that pD and pA are right

		//here have a list that gets a new entry with the tag indentification pair

		return newCode;
	}

	


	// DEBUG CENTER POINT
	if (debug == true) {
		circle(image, QRPos, 2, Scalar(0, 0, 255), -1, 8, 0);
		if (found) {
			try {
				float dist = cv_vectorSize(prevPos - QRPos);
				if (dist>=150.0)
					cout << dist << endl;
			}
			catch (const std::exception&) {}

			prevPos = QRPos;
		}
	}

	QRCode test(0);
	return test;
}

//Horrible Performance
int cv_HarrisCorner(Mat img) {

	cvtColor(img, img, CV_RGB2GRAY);
	cornerHarris(img, img, 2, 3, 0.05);
	imshow("Harris", img);

	return 0;
}
//###############################################################################
//Decodation Support Function
//###############################################################################
int decode(Mat inputImage) {
	ImageScanner scanner;
	scanner.set_config(ZBAR_NONE, ZBAR_CFG_ENABLE, 1);
	// obtain image data  
	//char file[256];
	//cin >> file;
	//Mat img = imread(file, 0);
	Mat imgout;
	//cvtColor(img, imgout, CV_GRAY2RGB);
	int width = inputImage.cols;
	int height = inputImage.rows;
	uchar *raw = (uchar *)inputImage.data;
	// wrap image data  
	Image imageFile(width, height, "Y800", raw, width * height);
	// scan the image for barcodes  
	int n = scanner.scan(imageFile);
	// extract results  
	for (Image::SymbolIterator symbol = imageFile.symbol_begin();
	symbol != imageFile.symbol_end();
		++symbol) {
		vector<Point> vp;
		// do something useful with results  
		cout << "decoded " << symbol->get_type_name()
			<< " symbol \"" << symbol->get_data() << '"' << " " << endl;
		int n = symbol->get_location_size();
		for (int i = 0; i<n; i++) {
			vp.push_back(Point(symbol->get_location_x(i), symbol->get_location_y(i)));
		}
		RotatedRect r = minAreaRect(vp);
		Point2f pts[4];
		r.points(pts);
		for (int i = 0; i<4; i++) {
			line(imgout, pts[i], pts[(i + 1) % 4], Scalar(255, 0, 0), 3);
		}
		cout << "Angle: " << r.angle << endl;
	}
	//imshow("imgout.jpg", imgout);

	//waitKey(0);

	return 1;
}


//###############################################################################
//Detection Support Functions 
//###############################################################################
int cv_findCorners(Point& pA, FiP fip_B, FiP fip_C, Point& pD, Point QRPos) {
	//gets the two Diagonal fips - bother other cornerPoints will be reconstructed from this
	//if FiP_A is know the outermost Point should be given in
	int decrementB = -1;
	int incrementB = 1;
	int decrementC = -1;
	int incrementC = 1;
	int indexOfPointB;
	Point cornerB = cv_getOuterCorner(fip_B, QRPos, indexOfPointB);
	int indexOfPointC;
	Point cornerC = cv_getOuterCorner(fip_C, QRPos, indexOfPointC);
	if (indexOfPointB == 0) {
		decrementB = 3;
	}
	else if (indexOfPointB == 3) {
		incrementB = -1;
	}
	if (indexOfPointC == 0) {
		decrementC = 3;
	}
	else if (indexOfPointC == 3) {
		incrementC = -1;
	}
	if (pA.x != 0 && pA.y != 0) { //This is stupid TODO check how properly check for null objects
		Point linePointB;
		Point linePointC;
		//the Point further away from the corner A will be the point with which a line for the Corner D can be built
		if (cv_euclideanDist(pA, fip_B.shape[indexOfPointB + decrementB]) > cv_euclideanDist(pA, fip_B.shape[indexOfPointB + incrementB]))
			linePointB = fip_B.shape[indexOfPointB + decrementB];
		else
			linePointB = fip_B.shape[indexOfPointB + incrementB];

		//the Point further away from the corner A will be the point with which a line for the Corner D can be built
		if (cv_euclideanDist(pA, fip_C.shape[indexOfPointC + decrementC]) > cv_euclideanDist(pA, fip_C.shape[indexOfPointC + incrementC]))
			linePointC = fip_C.shape[indexOfPointC + decrementC];
		else
			linePointC = fip_C.shape[indexOfPointC + incrementC];

		cv_getIntersection(cornerB, linePointB, cornerC, linePointC, pD);
		//cout << pD << endl;
		//vector<Point> lnB;
		//vector<Point> lnC;
		//lnB.push_back(cornerB);
		//lnB.push_back(linePointB);
		//lnC.push_back(cornerC);
		//lnC.push_back(linePointC);
		//drawContours(image, vector<vector<Point> >(1, lnB), -1, Scalar(255, 0, 0), 3, 8);
		//drawContours(image, vector<vector<Point> >(1, lnC), -1, Scalar(255, 255, 0), 3, 8);

	}
	else { // if only diagonal FiPs a know we have to cross both lines
		//generate all 4 lines
		Point linePointB1 = fip_B.shape[indexOfPointB + decrementB];
		Point linePointB2 = fip_B.shape[indexOfPointB + incrementB];
		Point linePointC1 = fip_C.shape[indexOfPointC + decrementC];
		Point linePointC2 = fip_C.shape[indexOfPointC + incrementC];

		//find angle to see which are roughly parallel and which are not
		float angle = cv_lineLineAngle(linePointB1, cornerB, linePointC1, cornerC); //Just try two if these are the parallel ones you know the two
																					 //Pairs that cross
		if (angle < 0.1) {
			//both are parallel
			cv_getIntersection(cornerB, linePointB2, cornerC, linePointC1, pD);
			cv_getIntersection(cornerB, linePointB1, cornerC, linePointC2, pA);
		}
		else {
			//both crossed
			cv_getIntersection(cornerB, linePointB1, cornerC, linePointC1, pD);
			cv_getIntersection(cornerB, linePointB2, cornerC, linePointC2, pA);
		}

		//try find two lines that cross for point (A)? how to check which point?
		
		//use other two to find the second point
	}

	return 0;
}

Point cv_getOuterCorner(FiP fip, Point center) {
	//gets the Point of the fip that makes up the Corner of the QR
	//should be Point furthest away from center
	Point outer(0, 0);
	float maxDist = 0.0;
	float dist;
	for (auto &p : fip.shape) {
		dist = cv_euclideanDist(p, center);
		if (dist >= maxDist) {
			maxDist = dist;
			outer = p;
		}
	}
	return outer;
}

Point cv_getOuterCorner(FiP fip, Point center, int& index) {
	//gets the Point of the fip that makes up the Corner of the QR
	//should be Point furthest away from center
	//this also returns the index of the Point in the FiP shape
	Point outer(0, 0);
	float maxDist = 0.0;
	float dist;
	for (int i = 0; i < 4; i++) {
		dist = cv_euclideanDist(fip.shape[i], center);
		if (dist >= maxDist) {
			maxDist = dist;
			outer = fip.shape[i];
			index = i;
		}
	}
	return outer;
}


vector<FiP> cv_getFiPOrder(vector<FiP> unordered){ //Returns the FiPs in order with 0 being A 1 being B 2 being C 
										 //if less then 3 FiP unknown FiPs should get the tag -1 so we can reconstruct them
	Point peak;
	Point hyotenuse1, hyotenuse2;

	// if three fips exist we can predict which is A B C easily by checking which Fips build the hypotenuse of the triangle
	if (unordered.size() == 3) {
		float side1 = cv_euclideanDist(unordered[0].pos, unordered[1].pos);
		float side2 = cv_euclideanDist(unordered[0].pos, unordered[2].pos);
		float side3 = cv_euclideanDist(unordered[1].pos, unordered[2].pos);

		float maxVal = max({ side1, side2, side3 });

		//This could be done better? Case statement maybe ?
		if (side1 == maxVal) {
			peak = unordered[2].pos;
			unordered[0].relPos = 1;
			unordered[1].relPos = 2;
			unordered[2].relPos = 0;
		}
		else if (side2 == maxVal) {
			peak = unordered[1].pos;
			unordered[0].relPos = 1;
			unordered[1].relPos = 0;
			unordered[2].relPos = 2;
		}
		else if (side3 == maxVal) {
			peak = unordered[0].pos;
			unordered[0].relPos = 0;
			unordered[1].relPos = 1;
			unordered[2].relPos = 2;
		}
	}
	// if two fips exist we can check if they are on the Diagonal and thus B and C or not and thus A and ?
	if (unordered.size() == 2) {
		//get the conecting side
		float side1 = cv_euclideanDist(unordered[0].pos, unordered[1].pos);
		//check if diagonal or orthogonal
		float angle = cv_lineLineAngle(unordered[0].shape[0], unordered[0].shape[1], unordered[0].pos, unordered[1].pos);
		if (0.1 < abs(angle) && abs(angle) < 0.9) {
			//diagonal
			//Point hypotenusePoint = (centers[0] + centers[1])*.5;
			unordered[0].relPos = 1;
			unordered[1].relPos = 2;
			//cout << "diagonal" << endl;
		}
		else {
			unordered[0].relPos = 0;
			unordered[1].relPos = 1;
			//cout << "parallel" << endl;
			//orthagonal/parallel
		}
	}

	// if one fip exist we can't build any order
	if (unordered.size() == 1) {
		unordered[0].relPos = 0;
	}

	return unordered;
}

vector<FiP> cv_reconstructMissingFiP(vector<FiP> orderedFiP) {
	
	vector<FiP> reconstructed;


	return reconstructed;
}

vector<Point> cv_getCornerPoints(vector<FiP> orderedReconstructedFiP) {
														//gets FiPs already in order and should return the points that make up the
														//Corners of the QR Codes
														//Returns a vector of Points making up the full Squares
	vector<Point> QRshape;


	return QRshape;
}


//###############################################################################
//General Support Functions 
//###############################################################################
bool cv_getIntersection(Point a1, Point a2, Point b1, Point b2, Point& res) {
	//Gets the intersectionPoint from 2 lines denoted by a1 a2 and b1 b2
	/*Point x = a2 - a1;
	Point d1 = b1 - a1;
	Point d2 = b2 - a2;

	float cross = (d1.x*d2.y) - (d1.y*d2.x);
	if (abs(cross) < /*EPS*//*1e-8) {
		cout << "FALSE" << endl;
		return false;
	}
		
	double t1 = ((x.x * d2.y) - (x.y * d2.x)) / cross;
	r = a1 + (d1 * t1);
	return true;*/
	Point p = a1;
	Point q = b1;
	Point r(a2 - a1);
	Point s(b2 - b1);

	if (cross(r, s) == 0) { return false; }

	float t = cross(q - p, s) / cross(r, s);

	res = p + t*r;
	return true;
}

float cross(Point2f v1, Point2f v2)
{
	return v1.x*v2.y - v1.y*v2.x;
}

bool cv_inRegion(Point center, int radius, Point newPoint) { //gives out true if in the region around the center is the newPoint
	if (cv_vectorSize(center- newPoint) <= radius) 
		return true;
	else
		return false;
}

float cv_lineLineAngle(Point l1_1, Point l1_2, Point l2_1, Point l2_2) { //Gets the angle that is between two lines denoted by their start and end Points
	//return ((l1_2.x - l1_1.x)*(l2_2.x - l2_1.x)) / (abs(l1_2.x - l1_1.x)*abs(l2_2.x - l2_1.x)); //<-- please test http://mathworld.wolfram.com/Line-LineAngle.html
	return ((l1_2 - l1_1).dot(l2_2 - l2_1)) / (cv_vectorSize(l1_2 - l1_1)*cv_vectorSize(l2_2 - l2_1));
}

float cv_vectorSize(Point a) {
	return hypot(a.x, a.y);
}


String cv_getOrientation(Point a, Point b) { //Gets orientation of Point b in relation to a
	String orientation;

	if (a.x >= b.x) {
		if (a.y >= b.y) {
			orientation = "se";
		}
		else {
			orientation = "ne";
		}
	}
	else {
		if (a.y >= b.y) {
			orientation = "sw";
		}
		else {
			orientation = "nw";
		}
	}

	return orientation;
}

Point cv_getCentroid(vector<Point> contour) {
	int sumX = 0, sumY = 0;
	int size = contour.size();
	Point centroid;
	if (size > 0) {

		for each(Point point in contour) {
			sumX += point.x;
			sumY += point.y;
		}

		centroid.x = sumX / size;
		centroid.y = sumY / size;
	}
	return centroid;
}

float cv_euclideanDist(Point p, Point q) {
	Point diff = p - q;
	return	hypot(diff.x, diff.y);
	//cv::sqrt(diff.x*diff.x + diff.y*diff.y);
}

bool cv_inRect(vector<Point> rectangle, Point p) { //ToDo: need knowledge in which order the Points are... do we have this?
	int minx, miny, maxx, maxy;
	if (rectangle.size() > 4) return false;
	/*for (int i = 0; i < 4; i++) { // show Corner Points
	circle(image, rectangle[i], 2, Scalar(255, 255, 0), -1, 8, 0);
	}*/
	minx = min({ rectangle[0].x, rectangle[1].x, rectangle[2].x, rectangle[3].x });
	miny = min({ rectangle[0].y, rectangle[1].y, rectangle[2].y, rectangle[3].y });
	maxx = max({ rectangle[0].x, rectangle[1].x, rectangle[2].x, rectangle[3].x });
	maxy = max({ rectangle[0].y, rectangle[1].y, rectangle[2].y, rectangle[3].y });
	if (p.x >= minx && p.x <= maxx && p.y >= miny && p.y <= maxy)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//For testing checks if Contour Center is in any Contour -- makes it easier to refresh  -- check which should be used later
bool cv_inFiPRegTesting(vector<vector<Point>>& FiPRegs, vector<Point> Contour, vector<bool>& updated) {
	Point p = cv_getCentroid(Contour);
	//circle(image, p, 2, Scalar(255, 0, 0), -1, 8, 0);
	int index = 0;
	if (FiPRegs.empty()) return false;
	for (std::vector<vector<Point>>::iterator it = FiPRegs.begin(); it != FiPRegs.end(); ++it) {
		//cout << "NUMBER :  " << index << "\n";
		if (cv_inRect(*it, p)) {
			// do the refresh here for testing later move it out but give back the index so we can do a different refresh for candidates
			//vector<Point> test = FiPRegs[index];// = Contour; // solve with pointers?
			FiPRegs[index] = Contour;
			updated[index].flip();// = true;
								  //*it = Contour;
								  //vector<Point> test2 = FiPRegs.at[index][0];
								  //cout << test[0];
								  //cout << test2[0];
			return true;
		}
		else {
			index++;
		}
	}
	return false;
}

//Check if ANY of the input contours center point is in the input Contour and returns the index of the first Point that fits // or the point? not sure yet
int cv_CandidateInRegion(vector<Point> contour, vector<vector<Point> > candidates) {
	int index = 0; // does it make sense to use a iterator when we need the index anyways? its also applies to the other two functions over this
	for (std::vector<vector<Point>>::iterator it = candidates.begin(); it != candidates.end(); ++it) {
		Point p = cv_getCentroid(*it);
		if (cv_inRect(contour, p)) {
			return index;
		}
	}
	return -1;
}

//Old Method that loaded all Pictures into Memory -- bit faster so maybe of use again at some point
/*
int imageLoadTest() {
vector<Mat> images;
Mat image;
String number = "";
bool moreFiles = true;
int i = 0;
int j = 0;
int key = 0;
vector<int> timePassedAll;
vector<float> fpsAll;
// Creation of Intermediate 'Image' Objects required later
Mat gray(image.size(), CV_MAKETYPE(image.depth(), 1));			// To hold Grayscale Image
Mat edges(image.size(), CV_MAKETYPE(image.depth(), 1));			// To hold Grayscale Image
Mat traces(image.size(), CV_8UC3);								// For Debug Visuals
Mat empty(Size(100,100), CV_MAKETYPE(image.depth(), 1));

vector<vector<Point> > contours;
vector<vector<Point> > fiPReg;
vector<Vec4i> hierarchy;
vector<Point> pointsseq;    //used to save the approximated sides of each contour

while (moreFiles) {
number = to_string(i);
while (number.size() < 7) {
number = "0" + number;
}
cout << number << endl;
image = imread("C:/Users/Frederik/Desktop/VidTests/720-frame-png/image-" + number + ".png");

if (image.data == NULL) {
cerr << "No More Files" << endl;
moreFiles = false;
}
else {
images.push_back(image);
}
i++;
}

chronoNow = chrono::system_clock::now();

while (key != 'q') //Main Computation Loop
{

chronoPrev = chronoNow;
chronoNow = chrono::system_clock::now();

image = images.at(j);


//Image Processing
cvtColor(image, gray, CV_RGB2GRAY);		// Convert Image captured from Image Input to GrayScale
//medianBlur(gray, gray, 3);
Canny(gray, edges, 150, 200, 3);		// Apply Canny edge detection on the gray image
findContours(edges, contours, hierarchy, RETR_TREE, CHAIN_APPROX_TC89_KCOS);

//Time Calculation
int timePassed = std::chrono::duration_cast<std::chrono::milliseconds>(chronoNow - chronoPrev).count();
float fps = 1000 / (float)timePassed;
timePassedAll.push_back(timePassed);
fpsAll.push_back(fps);
//cout << "Time passed since last frame" << timePassed << "\n";
//cout << "Frames per second " << fps << "\n";

if (showCalc)
imshow("image", image);
else
imshow("image", empty);

j++;
if (j == (i - 1)) {
break;
}
key = waitKey(1);
}

//Calculate FPS and give a history of time taken for Webcam
float fpsSum = 0;
for (int i = 0; i < timePassedAll.size(); i++) {
//cout << timePassedAll.at(i) << " - " ;
if (i != 0) {
fpsSum = fpsSum + fpsAll.at(i);
}
}
//cout << "\n";
float fpsAverage = fpsSum / (float)(fpsAll.size() - 1);
cout << "Average FPS: " << fpsAverage << "\n";

QRTest();
return -1;
}
*/
