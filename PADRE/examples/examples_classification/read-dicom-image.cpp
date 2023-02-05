#include <iostream>
#include "dicom/DicomReader.h"
#include "svm/svm.h"

using namespace std;

// Read a dicom image and print its image array

int main()
{
  DicomReader dicomObj("img/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");
  int **arr = dicomObj.getImageArray(12);
  for (int i = 0; i < dicomObj.getHeight(); ++i)
  {
    for (int j = 0; j < dicomObj.getWidth(); ++j)
    {
      cout << arr[i][j] << " ";
    }
  }

  return 0;
}
