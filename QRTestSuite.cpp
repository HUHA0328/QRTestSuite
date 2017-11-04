// QRTestSuite.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "QRTestSuite.hpp"

//###############################################################################
//Methods
//###############################################################################
QRFinder::QRFinder(int inputCameraImageWidth,
	int inputCameraImageHeight,
	const cv::Mat_<double>& inputCameraCalibrationMatrix,
	const cv::Mat_<double>& inputCameraDistortionParameters,
	bool pShowResults) 
{


	// Check inputs
	if (inputCameraImageWidth <= 0 || inputCameraImageHeight <= 0) {
		std::cerr << "Camera image dimensions invalid!" << std::endl;
	}

	// Check camera calibration matrix dimensions
	if (inputCameraCalibrationMatrix.dims != 2) {
		std::cerr << "Camera calibration matrix is not 3x3!" << std::endl;
	}

	for (int i = 0; i < inputCameraCalibrationMatrix.dims; i++) {
		if (inputCameraCalibrationMatrix.size[i] != 3) {
			std::cerr << "Camera calibration matrix is not 3x3!" << std::endl;
		}
	}

	// Check distortion coefficients size
	if (inputCameraDistortionParameters.dims != 2) {
		std::cerr << "Distortion coefficents vector is not 5x1!" << std::endl;
	}

	if (inputCameraDistortionParameters.size[1] != 1 ||
		inputCameraDistortionParameters.size[0] != 5) {
		std::cerr << "Distortion coefficents vector is not 5x1!" << std::endl;
	}

	expectedCameraImageWidth = inputCameraImageWidth;
	expectedCameraImageHeight = inputCameraImageHeight;

	cameraMatrix = inputCameraCalibrationMatrix;
	distortionParameters = inputCameraDistortionParameters;
	showResults = pShowResults;

	// Configure the QR code reader object
	zbarScanner.set_config(zbar::ZBAR_QRCODE, zbar::ZBAR_CFG_ENABLE, 1);
	zbarScanner.enable_cache(false);  // Set it so that it will show QR code
									  // result even if it was in the last frame

									  // Create window to show results, if we are suppose to
	if (showResults == true) {
		cv::namedWindow(QRCodeStateEstimatorWindowTitle, CV_WINDOW_AUTOSIZE);
	}
};



class QRFinder::FiP {
public:
	FiP ();
	FiP (int pRelPos);
	FiP (Point pPos, vector<Point> pShape);
	Point pos;
	vector<Point> shape;
	int relPos; //relative pos in the QR code: 0 the corner, 1 and 2 the diagonal Fips: unkown -1? or 0?

	void setRelPos(int pRelPos) {
		relPos = pRelPos;
	}
};

QRFinder::FiP::FiP() {
	
}

QRFinder::FiP::FiP(int pRelPos) {
	relPos = pRelPos;
}

QRFinder::FiP::FiP (Point pPos, vector<Point> pShape) {
	pos = pPos;
	shape = pShape;
}

class QRFinder::QRCode { //Class for the QR-Code that saves all relevant information for one Instance
	public:
		QRCode(int ptag);
		Point pos; // x and y pos of the middle of the QR-Code (mostly for testing as the Corners will be used for processing)
		int tag; // Tag that will identify this QR-Code so that every QR-data will only have one corresponding QR-Code
		FiP fip_a, fip_b, fip_c; //The finderpatterns
		Point a, b, c, d; //the Corner Points
		String orientation; //does it make sense to create a new class for this?
		bool decode_success; //shows it was successfull decoded
		String data; //gives the data -> maybe we want to sace this in a table with the tags so a possible parallel process wouldn't have to access the Object everytime
		//int screenSize; //Size in ScreenSpace
		//float realSize; //Size in Centimeters - probably 20cm 

};

QRFinder::QRCode::QRCode(int ptag) {
	tag = ptag;
}

int QRFinder::QRBenchmark(bool showCalc)
{

	int64 e1, e2;
	float t;
	float globalFPS = 0.0;
	vector<int> timePassedAll;
	vector<vector<FiP>> FiPList;
	vector<QRCode> QRList;

	//capture.open("C:/Users/Frederik/Desktop/VidTests/long-exposure.mp4");
	//capture.open("C:/Users/Frederik/Desktop/VidTests/720p-dist-move.mp4");
	//capture.open("C:/Users/Frederik/Desktop/VidTests/long-exposure-HQ.mp4");
	//capture.open("C:/Users/Frederik/Desktop/VidTests/moto/single-short-leaving.mp4");
	//capture.open("C:/Users/Frederik/Desktop/VidTests/moto/single-short-easy.mp4");
	//capture.open("C:/Users/Frederik/Desktop/VidTests/moto/single-long-all.mp4");
	//capture.open("C:/Users/Frederik/Desktop/VidTests/moto/multi-short-distance.mp4");
	VideoCapture capture;
	capture.open("C:/Users/Frederik/Desktop/VidTests/moto/multi-short-movement.mp4");

	// Creation of Intermediate 'Image' Objects required later
	Mat empty(Size(100, 100), CV_MAKETYPE(image.depth(), 1));

	if (!capture.isOpened()) {
		cerr << " ERR: Unable find input Video source." << endl;
		return -1;
	}

	Mat empty(Size(100, 100), CV_MAKETYPE(image.depth(), 1));
	e1 = getTickCount();

	while (capture.isOpened()) //Main Computation Loop
	{

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
			/*cout << "-----------Accuracy Metrics-----------" << endl;
			cout << "Frames solo           : " << ((float)solo / (float)frames) * 100 << "%" << endl;  //one found
			cout << "Frames partial        : " << ((float)partial / (float)frames) * 100 << "%" << endl;  //two found
			cout << "Frames full detection : " << ((float)full / (float)frames) * 100 << "%" << endl; //Three FiPs corrected or Reconstructed
			cout << "Frames overshoot      : " << over << endl;
			cout << "--------------QR decoded--------------" << endl;
			cout << "QRCode identified     : " << ((float)decoded / (float)frames) * 100 << "%" << endl;
			cout << "tags given out        : " << tagDataMap.size() << endl;*/

			return 0;
		}

		//Image Processing
		FiPList = cv_FiPdetection(image, FiPList); //<--- the QR list should go there

												   //cout << FiPList.size() << endl;

												   //grouping 

		for (auto &fip : FiPList) {
			QRCode qrNew = cv_QRdetection(fip, QRList);
			if (qrNew.tag != 0) {
				for (int i = 0; i < QRList.size(); i++) { //<- It should also be possible to optimize this but since we technicly should never get more than exsisting QR codes (so around 3) entrys it should always run very fast
					if (QRList[i].tag == qrNew.tag) {
						QRList.erase(QRList.begin() + i);
						break;
					}
				}
				QRList.push_back(qrNew);
			}
		}
		
		//Debug Rendering
		if (showCalc)
			imshow("image", image);
		else
			imshow("image", empty);

		waitKey(1);
	}
}

int QRFinder::QRScan(vector<vector<FiP>> FiPList, vector<QRCode> QRList,
	const cv::Mat& inputGrayscaleFrame,
	cv::Mat& inputCameraPoseBuffer,
	std::string& inputQRCodeIdentifierBuffer,
	double& inputQRCodeDimensionBuffer)
{
	FiPList = cv_FiPdetection(image, FiPList); //<--- the QR list should go there

											   //cout << FiPList.size() << endl;

											   //grouping 

	for (auto &fip : FiPList) {
		QRCode qrNew = cv_QRdetection(fip, QRList);
		if (qrNew.tag != 0) {
			for (int i = 0; i < QRList.size(); i++) { //<- It should also be possible to optimize this but since we technicly should never get more than exsisting QR codes (so around 3) entrys it should always run very fast
				if (QRList[i].tag == qrNew.tag) {
					QRList.erase(QRList.begin() + i);
					break;
				}
			}
			QRList.push_back(qrNew);
		}
	}
	
}

/* calling a method an reinitializing a lot of Variables should create a overhead compared to just use one method?
	Maybe I should do a speed test since it would be far easier to understand in more methods but if it costs more time
	nothing is really gained*/
vector<vector<QRFinder::FiP>> QRFinder::cv_FiPdetection(Mat inputImage, vector<vector<FiP>> prevImage ) /* //<--- the QR list should go there
									input:

									output:

									*/
{

	vector<vector<FiP>> fipList;
	vector<vector<Point> > contours, contoursMorph;
	vector<Vec4i> hierarchy, hierarchyMorph;
	vector<Point> pointsseq;    //used to save the approximated sides of each contour
	vector<vector<Point> > finderPatterns;
	vector<vector<Point> > finderCandidate;
	vector<bool> updated;
	Mat morphImage;
	updated.clear();

	

	//cout << "test" << endl;
	cvtColor(inputImage, inputImage, CV_RGB2GRAY);										// Convert Image captured from Image Input to GrayScale
	resize(inputImage, inputImage, Size(inputImage.size().width / 2, inputImage.size().height / 2));
	//resize(inputImage, inputImage, Size(640, 360));
	threshold(inputImage, inputImage, 180, 255, THRESH_BINARY);							//<- Probably some kind of local threshold better in the Areas of Interest 
	//imshow("debug", inputImage);
	morphologyEx(inputImage, morphImage, MORPH_OPEN, getStructuringElement(MORPH_RECT, Size(6, 6)), Point(-1,-1), 3); // these behave oppsoite to what to expect since white is obviosly value even if for this purpose it signals that there is nothing
	morphologyEx(morphImage, morphImage, MORPH_CLOSE, getStructuringElement(MORPH_RECT, Size(5, 5)), Point(-1, -1), 3);
	resize(inputImage, inputImage, Size(inputImage.size().width * 2, inputImage.size().height * 2));
	resize(morphImage, morphImage, Size(morphImage.size().width * 2, morphImage.size().height * 2));
	copyMakeBorder(morphImage, morphImage, 10, 10, 10, 10, BORDER_CONSTANT, Scalar(255, 255, 255)); //<---- YES! But adjust for padding
	copyMakeBorder(inputImage, inputImage, 10, 10, 10, 10, BORDER_CONSTANT, Scalar(255, 255, 255)); //<---- for testing
	copyMakeBorder(image, image, 10, 10, 10, 10, BORDER_CONSTANT, Scalar(255, 255, 255)); //<---- for testing

	//imshow("morph", morphImage);
	//resize(inputImage, inputImage, Size(1920, 1080));
	//waitKey(0);
	//imshow("debug", inputImage);
	//imshow("morph", morphImage);



	
	//Canny(inputImage, inputImage, 80, 150, 3);											// Apply Canny edge detection on the gray image
	//findContours(inputImage, contours, hierarchy, RETR_TREE, CHAIN_APPROX_TC89_KCOS);
	Canny(morphImage, morphImage, 80, 150, 3);											// Apply Canny edge detection on the gray image
	findContours(morphImage, contoursMorph, hierarchyMorph, RETR_TREE, CHAIN_APPROX_TC89_KCOS);
	//imshow("debugContours", morphImage);
	//cout << "before " << contours.size() << endl;
	//cout << "after " << contoursMorph.size() << endl;
	//cout << ((float)contoursMorph.size() / (float)contours.size()) * 100 << "%" << endl;
	//cout << endl;
	vector<Mat> image_roi;
	int j = 0;
	for (int i1 = 0; i1 < contoursMorph.size(); i1 += 2) {
		approxPolyDP(contoursMorph[i1], pointsseq, arcLength(contoursMorph[i1], true)*0.05, true);
		if (pointsseq.size() == 4 && isContourConvex(pointsseq)) {
			vector<FiP> qrCandidate;
			vector<vector<Point> > fiPReg;
			//Maybe jump over obviously no QR-Codes ?
			Rect regionOfInterest = cv_getRect(pointsseq[0], pointsseq[1], pointsseq[2], pointsseq[3], true);
			image_roi.push_back(image(regionOfInterest));
			cvtColor(image_roi[j], image_roi[j], CV_RGB2GRAY);
			threshold(image_roi[j], image_roi[j], 180, 255, THRESH_BINARY);
			copyMakeBorder(image_roi[j], image_roi[j], 5, 5, 5, 5, BORDER_CONSTANT, Scalar(255, 255, 255));
			resize(image_roi[j], image_roi[j], Size(image_roi[j].size().width * 2, image_roi[j].size().height * 2));
			Canny(image_roi[j], image_roi[j], 10, 200, 3);// 80, 150, 3);
			findContours(image_roi[j], contours, hierarchy, RETR_TREE, CHAIN_APPROX_TC89_KCOS);
			j++;

			for (int i2 = 0; i2 < contours.size(); i2++)
			{
				//Find the approximated polygon of the contour we are examining
				approxPolyDP(contours[i2], pointsseq, arcLength(contours[i2], true)*0.02, true); // <- this somehow doesnt work
				if (pointsseq.size() == 4)      // only quadrilaterals contours are examined // <- this neither or I misunderstand it
				{
					int k = i2;
					int c = 0;

					while (hierarchy[k][2] != -1)
					{
						k = hierarchy[k][2];
						c = c + 1;
					}

					if (c >= 4)
					{
						finderPatterns.push_back(contours[k]);
						vector<Point> fipSquare;
						approxPolyDP(contours[k - c], fipSquare, arcLength(contours[i2], true)*0.02, true);
						//fipSquare = contours[k - c];
						for (int i3 = 0; i3 < fipSquare.size(); i3++) {
							fipSquare[i3].x = fipSquare[i3].x/2 + regionOfInterest.x -5;
							fipSquare[i3].y = fipSquare[i3].y/2 + regionOfInterest.y -5;
						}
						fiPReg.push_back(fipSquare);
						if (c >= 5) {
							i2 += 2;
						}
					} 
					else
						i2 += c;
				}
			}
			for (auto &fip : fiPReg) {
				Point Center = cv_getCentroid(fip);
				qrCandidate.push_back(FiP(Center, fip));
			}
			if (qrCandidate.size()>0)
				fipList.push_back(qrCandidate);
		}

	}
	

	/*
	for (int z = 0; z < fipList.size(); z++) {
		Scalar groupCol;
		if (z == 0)
			groupCol = Scalar(255, 0, 0);
		else if (z == 1) {
			groupCol = Scalar(0, 255, 0);
			//waitKey(0);
		}
		else if (z == 2)
			groupCol = Scalar(255, 255, 0);
		else
			groupCol = Scalar(0, 0, 255);
		for (int zz = 0; zz < fipList[z].size(); zz++) {
			drawContours(image, vector<vector<Point> >(1, fipList[z][zz].shape), -1, groupCol, 5, 8);
		}
	}*/
	
	

	//drawContours(image, fiPReg, -1, Scalar(0, 0, 255), 1, 8);
	return fipList;
}

QRFinder::QRCode QRFinder::cv_QRdetection(vector<FiP> fipImage, vector<QRCode> qrPrevImage) {
	QRCode returnCode(0);
	Point QRPos;
	FiP fip_A, fip_B, fip_C;
	bool found = false;
	bool success = false;
	String qrData = "";
	int tag;
	int oldQRindex;

	// process FiPs
	fipImage = cv_getFiPOrder(fipImage);
	// if 3 or 2 diagonal we already have the position
	//cout << fipImage.size() << endl;

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
			fip_A = fipImage[2]; //2 //vorher 1
			fip_B = fipImage[1]; //1 //2
			fip_C = fipImage[0]; //0 //0
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
		else {
			//Two parallel FiPs
			fip_A = fipImage[0];
			fip_B = fipImage[1];
			QRCode oldQR = QRCode(0);
			float mindist = 9999;
			int i = 0;
			for (auto &qr : qrPrevImage) {
				float newDist = cv_euclideanDist(fip_A.pos, qr.pos);
				if (newDist < mindist) {
					QRPos = qr.pos;
					mindist = newDist;
					oldQR = qr;
					oldQRindex = i;
				}
				i++;
			}
			QRPos = cv_find2FipPos(fip_A, fip_B, oldQR); // this doesn't work currently since it will not exclude other QRs in the image right now || it Should work once the code was properly found and there wasn't a big movement in between
			//SO IT WILL BE OVERWRITTEN BY THE NOT FOUND CLAUSE 

			//found = true;
		}
	}
	else if (fipImage.size() == 1) {
		// if 1 check possible positions
		// What could we do?
		// check if there is a tag expected in the radius?
		// ignore might work best? (only if this happens rarely which seems to be the case)
		fip_A = fipImage[0];
		//QRPos = cv_find1FiPPos(fip_A, qrPrevImage);
		float mindist = 9999;
		for (auto &qr : qrPrevImage) {
			float newDist = cv_euclideanDist(fip_A.pos, qr.pos);
			if (newDist < mindist) {
				QRPos = qr.pos;
				mindist = newDist;
			}
		}
		//QRPos = qrPrevImage.pos;
	}

	// check if tag exists else give new one
	// we want a distance to a tag from a previous image to check if this is the same tag
	// this could create problems when there are overlaying tags
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
		//cout << "FiPs found" << fipImage.size() << endl;
		//qrPrevImage.pos = QRPos;
		//If none is found return a QRCode with tag 0
		returnCode = QRCode(0);
	}
	else {
		QRCode oldQR = QRCode(0);
		float mindist = 9999;
		for (auto &qr : qrPrevImage) {
			float newDist = cv_euclideanDist(QRPos, qr.pos);
			if (newDist < mindist) {
				mindist = newDist;
				oldQR = qr;
			}
		}
		if (mindist <= 150.0 && oldQR.decode_success == true) {
			// check if there is tag for the Pos
			// if yes load data from app

			//if (cv_vectorSize(qrPrevImage.pos - QRPos) <= 150.0) {
			//cout << "direct" << endl;
			oldQR.pos = QRPos;
			returnCode = oldQR;
		}
		else {
			// if no reconstruct planar QR code and load data
			//Identify

			//get CornerPointA in Case it exists already
			Point pA;
			if (fip_A.relPos == -1)
				pA = Point(0, 0);
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
			qrImg.push_back(Point2f(cv_getOuterCorner(fip_B, QRPos)));
			qrImg.push_back(Point2f(pD)); //<- order important?
			qrImg.push_back(Point2f(cv_getOuterCorner(fip_C, QRPos)));

			//for drawining debug version
			vector<Point> qrImg2;
			qrImg2.push_back(Point(pA));
			qrImg2.push_back(Point(cv_getOuterCorner(fip_B, QRPos)));
			qrImg2.push_back(Point(pD)); //<- order important?
			qrImg2.push_back(Point(cv_getOuterCorner(fip_C, QRPos)));

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

			if (qrImg.size() == 4 && dst.size() == 4)			// Failsafe for WarpMatrix Calculation to have only 4 Points with src and dst
			{
				warp_matrix = getPerspectiveTransform(qrImg, dst);
				//warpPerspective(image, qr_raw, warp_matrix, Size(200, 200));
				warpPerspective(image, qr_raw, warp_matrix, Size(200, 200)); // we could get this from the already worked on image!
				copyMakeBorder(qr_raw, qr, 10, 10, 10, 10, BORDER_CONSTANT, Scalar(255, 255, 255));

				cvtColor(qr, qr_gray, CV_RGB2GRAY); //could save this

				//check image quality before processing?
				threshold(qr_gray, qr_thres, 180, 255, CV_THRESH_BINARY); //and this

				qrData = decode(qr_thres);
				
				if (qrData != "ERROR")
					success = true;
				else {
					success = false;
				}

			}

			if (success) {
				//If tag exists retag
				tag = findTaginList(qrData);
				if (tag == tagDataMap.size())
					tagDataMap.push_back(qrData);
				//int tag = tagnumber++; //<- usually find next free tag - global variable?

				//create QRCode Object
				QRCode newCode(tag);
				newCode.pos = QRPos;
				newCode.decode_success = success;
				newCode.data = qrData; // <-- we need this or always over tag?

				//if (fipImage.size() == 3) //<- only do this if we had all threeFiPs so its sure that pD and pA are right - not sure if this is the best idea anymore.
				newCode.orientation = cv_getOrientation(pA, pD);

				//here have a list that gets a new entry with the tag indentification pair

				returnCode = newCode;
			}
			else {
				// TODO : WHAT SHOULD HAPPEN HERE IN MULTI
				QRCode newCode(0);
				newCode.pos = QRPos;
				returnCode = newCode; // <-- will get deleted anywaysbecause tag == 0 
			}
		}
	}


	//cout << "THIS WILL NEVER BE REACHED" << endl;
	// DEBUG CENTER POINT
	if (showResults == true) {
		circle(image, QRPos, 10, Scalar(0, 0, 255), -1, 8, 0);
	}

	return returnCode;
}


//###############################################################################
//Decodation Support Function
//###############################################################################
String QRFinder::decode(Mat inputImage) {
	ImageScanner scanner;
	scanner.set_config(ZBAR_NONE, ZBAR_CFG_ENABLE, 1);
	bool found = false;
	String output;
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
		output = symbol->get_data();
		found = true;
	}

	if (found)
		return output;
	else
		return "ERROR";
}

String QRFinder::decodeMobile(Mat inputImage) {
	bool found = false;
	String output;
	Mat imgout;

	int width = inputImage.cols;
	int height = inputImage.rows;
	uchar *raw = (uchar *)inputImage.data;
	// wrap image data  
	Image imageFile(width, height, "Y800", raw, width * height);
	// scan the image for barcodes  
	int n = zbarScanner.scan(imageFile);
	// extract results  
	for (Image::SymbolIterator symbol = imageFile.symbol_begin();
	symbol != imageFile.symbol_end();
		++symbol) {
		output = symbol->get_data();
		found = true;
	}

	if (found)
		return output;
	else
		return "ERROR";
}


//###############################################################################
//Detection Support Functions 
//###############################################################################
Point QRFinder::cv_find1FiPPos(QRFinder::FiP fip_A, QRFinder::QRCode qrPrevImage) {
	String orientation = qrPrevImage.orientation;
	Point A1 = fip_A.shape[0];
	Point A2 = fip_A.shape[1];
	Point A3 = fip_A.shape[2];
	Point A4 = fip_A.shape[3];
	//Point distenceVectorPrev = qrPrevImage



	//if (orientation == cv_getOrientation(fip_A.shape[0], fip_A.shape[2]))
	return Point(0, 0); // <-debug
}

Point QRFinder::cv_find2FipPos(FiP fip_A, FiP fip_B, QRCode prevQR) {
	Point p_1;
	Point p_2;
	bool switched = false;
		
	Point pA_1 = fip_A.shape[0];
	Point pA_2 = fip_A.shape[2];

	Point pB_1 = fip_B.shape[1];
	Point pB_2 = fip_B.shape[3];

	//shouldnt happen BUT check if parallel if yes switch one
	float angle = cv_lineLineAngle(pA_1, pA_2, pB_1, pB_2); 
	
	
	if (angle > 0.1) {
		//both are parallel
		pB_1 = fip_B.shape[0];
		pB_2 = fip_B.shape[2];
		switched = true;
	}

	//they should Meet -- switch Both for other Point 
	cv_getIntersection(pA_1, pA_2, pB_1, pB_2, p_1);

	if (!switched) {
		pA_1 = fip_A.shape[1];
		pA_2 = fip_A.shape[3];

		pB_1 = fip_B.shape[0];
		pB_2 = fip_B.shape[2];
	}
	else {
		pA_1 = fip_A.shape[1];
		pA_2 = fip_A.shape[3];

		pB_1 = fip_B.shape[1];
		pB_2 = fip_B.shape[3];
	}

	cv_getIntersection(pA_1, pA_2, pB_1, pB_2, p_2);

	//get orientation from prevQR
	String orientation = prevQR.orientation; 

	//the Point that has the same orientation to one of the known FiPs is the right Point
	if (orientation == cv_getOrientation(pA_1, p_1))
		return p_1;
	else if (orientation == cv_getOrientation(pB_2, p_1))
		return p_1;
	else
		return p_2;
}

int QRFinder::cv_findCorners(Point& pA, FiP fip_B, FiP fip_C, Point& pD, Point QRPos) {
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
		incrementB = -3;
	}
	if (indexOfPointC == 0) {
		decrementC = 3;
	}
	else if (indexOfPointC == 3) {
		incrementC = -3;
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

Point QRFinder::cv_getOuterCorner(FiP fip, Point center) {
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

Point QRFinder::cv_getOuterCorner(FiP fip, Point center, int& index) {
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


vector<QRFinder::FiP> QRFinder::cv_getFiPOrder(vector<FiP> unordered){ //Returns the FiPs in order with 0 being A 1 being B 2 being C 
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
		//get the connecting side
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

//###############################################################################
//General Support Functions 
//###############################################################################
Rect QRFinder::cv_getRect(Point p1, Point p2, Point p3, Point p4, bool pad = false) {
	int minX, minY, maxX, maxY;
	
	minX = min({ p1.x, p2.x, p3.x, p4.x });
	maxX = max({ p1.x, p2.x, p3.x, p4.x });
	minY = min({ p1.y, p2.y, p3.y, p4.y });
	maxY = max({ p1.y, p2.y, p3.y, p4.y });
	if (pad == true) {
		minX = (int)minX*0.95;
		maxX = (int)maxX*1.05;
		minY = (int)minY*0.95;
		maxY = (int)maxY*1.05;
	}

	Rect newRect(Point(minX, minY), Point(min(maxX, image.size().width), min(maxY, image.size().height)));
	return newRect;
}

int QRFinder::findTaginList(String inputData) { //finds the tag if the QR code was already saved else gives out the next free tag 
	int tag = 0;

	for (int i = 0; i < tagDataMap.size(); i++) {
		if (tagDataMap[i] == inputData)
			return i;
	}
	return tagDataMap.size();
}

bool QRFinder::cv_getIntersection(Point a1, Point a2, Point b1, Point b2, Point& res) {
	//Gets the intersectionPoint from 2 lines denoted by a1 a2 and b1 b2
	Point p = a1;
	Point q = b1;
	Point r(a2 - a1);
	Point s(b2 - b1);

	if (cross(r, s) == 0) { return false; }

	float t = cross(q - p, s) / cross(r, s);

	res = p + t*r;
	return true;
}

float QRFinder::cross(Point2f v1, Point2f v2)
{
	return v1.x*v2.y - v1.y*v2.x;
}

bool QRFinder::cv_inRegion(Point center, int radius, Point newPoint) { //gives out true if in the region around the center is the newPoint
	if (cv_vectorSize(center- newPoint) <= radius) 
		return true;
	else
		return false;
}

float QRFinder::cv_lineLineAngle(Point l1_1, Point l1_2, Point l2_1, Point l2_2) { //Gets the angle that is between two lines denoted by their start and end Points
	// http://mathworld.wolfram.com/Line-LineAngle.html
	return ((l1_2 - l1_1).dot(l2_2 - l2_1)) / (cv_vectorSize(l1_2 - l1_1)*cv_vectorSize(l2_2 - l2_1));
}

float QRFinder::cv_vectorSize(Point a) {
	return hypot(a.x, a.y);
}


String QRFinder::cv_getOrientation(Point a, Point b) { //Gets orientation of Point b in relation to a
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

Point QRFinder::cv_getCentroid(vector<Point> contour) {
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

float QRFinder::cv_euclideanDist(Point p, Point q) {
	Point diff = p - q;
	return	hypot(diff.x, diff.y);
	//cv::sqrt(diff.x*diff.x + diff.y*diff.y);
}

bool QRFinder::cv_inRect(vector<Point> rectangle, Point p) { //ToDo: need knowledge in which order the Points are... do we have this?
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
bool QRFinder::cv_inFiPRegTesting(vector<vector<Point>>& FiPRegs, vector<Point> Contour, vector<bool>& updated) {
	Point p = cv_getCentroid(Contour);
	//circle(image, p, 10, Scalar(255, 0, 0), -1, 8, 0);
	int index = 0;
	if (FiPRegs.empty()) return false;
	for (std::vector<vector<Point>>::iterator it = FiPRegs.begin(); it != FiPRegs.end(); ++it) {

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
int QRFinder::cv_CandidateInRegion(vector<Point> contour, vector<vector<Point> > candidates) {
	int index = 0; // does it make sense to use a iterator when we need the index anyways? its also applies to the other two functions over this
	for (std::vector<vector<Point>>::iterator it = candidates.begin(); it != candidates.end(); ++it) {
		Point p = cv_getCentroid(*it);
		if (cv_inRect(contour, p)) {
			return index;
		}
	}
	return -1;
}