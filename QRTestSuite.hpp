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
int64 e1, e2;
float t;
vector<float> allFPS;
float avgFPS = 0.0, highFPS = 0.0, lowFPS = 999.0, globalFPS = 0.0; //FPS Performance Calculation happening global so minimal impact on algorithm
int decoded = 0;
vector<int> fips; //just saves the number of FiPs detected for Benchmark
Point prevPos;
//int tagnumber = 0;
vector<string> tagDataMap;

//#################		Images
Mat image; //back to multiFIPdetector

		   //#################		Classes
class FiP;
class QRCode;

//#################		Main Methods
int QRTest();
int fpsCalc(float fps);
//int imageLoadTest();

//Detector Methods
vector<vector<FiP>> cv_FiPdetection(Mat inputImage, vector<vector<FiP>> prevImage);
QRCode cv_QRdetection(vector<FiP> fipImage, vector<QRCode> qrpPrevImage);
int cv_HarrisCorner(Mat img);
vector<FiP> cv_getFiPOrder(vector<FiP> unordered);

//#################		Support Methods
Rect cv_getRect(Point p1, Point p2, Point p3, Point p4, bool pad);
int findTaginList(String inputData);
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
Point cv_find2FipPos(FiP fip_A, FiP fip_B, QRCode prevQR);
//ZBAR Decode Function
String decode(Mat inputImage);