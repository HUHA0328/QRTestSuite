#ifndef QRFINDERHPP
#define QRFINDERHPP

#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/videoio.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <chrono>
#include <zbar.h>
using namespace std;
using namespace cv;
using namespace zbar;
static const std::string QRCodeStateEstimatorWindowTitle =
"QR Code Finder";
//###############################################################################
// Initialization
//###############################################################################

class QRFinder {

	class FiP;
	class QRCode;

public:



	QRFinder(int inputCameraImageWidth,
		int inputCameraImageHeight,
		const cv::Mat_<double>& inputCameraCalibrationMatrix,
		const cv::Mat_<double>& inputCameraDistortionParameters,
		bool pShowResults) {};

	//#################		Global Variables
	vector<int> fips; //just saves the number of FiPs detected for Benchmark
	//Point prevPos;
	//int tagnumber = 0;
	vector<string> tagDataMap;
	int expectedCameraImageWidth;
	int expectedCameraImageHeight;
	cv::Mat_<double> cameraMatrix;  // 3x3 matrix
	cv::Mat_<double> distortionParameters;  // 1x5 matrix
	zbar::ImageScanner zbarScanner;
	cv::Mat frameBuffer;
	bool showResults; // True if the image should be shown in a window


	//#################		Images
	Mat image; //back to multiFIPdetector

	//#################		Classes
	//class FiP;
	//class QRCode;

	//#################		Main Methods
	int QRBenchmark(bool showCalc);
	int QRScan(
		vector<vector<FiP>> FiPList,
		vector<QRCode> QRList,
		const cv::Mat& inputGrayscaleFrame,
		cv::Mat& inputCameraPoseBuffer,
		std::string& inputQRCodeIdentifierBuffer,
		double& inputQRCodeDimensionBuffe);

	//Detector Methods
	vector<vector<FiP>> cv_FiPdetection(Mat inputImage, vector<vector<FiP>> prevImage);
	QRCode cv_QRdetection(vector<FiP> fipImage, vector<QRCode> qrpPrevImage);
	vector<FiP> cv_getFiPOrder(vector<FiP> unordered);
	int CalculatePoseFromQRList(const cv::Mat& inputFrame,
		vector<QRCode> QRList,
		std::vector<cv::Mat>& inputCameraPosesBuffer,
		std::vector<std::string>& inputQRCodeIdentifiersBuffer,
		std::vector<double>& inputQRCodeDimensionsBuffer);
	float getQRSizefromString(String qrString);

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
	float cross(Point2f v1, Point2f v2);
	String cv_getOrientation(Point a, Point b);
	Point cv_find2FipPos(FiP fip_A, FiP fip_B, QRCode prevQR);
	Point cv_find1FiPPos(FiP fip_A, QRCode qrPrevImage);
	bool cv_inRegion(Point center, int radius, Point newPoint);
	//ZBAR Decode Function
	String decode(Mat inputImage);
	String decodeMobile(Mat inputImage);
};

#endif
