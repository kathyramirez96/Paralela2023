#include "tiny_dnn/tiny_dnn.h"
#include "dicom/DicomReader.h"

using namespace tiny_dnn;
using namespace tiny_dnn::activation;
using namespace tiny_dnn::layers;
using namespace std;

//    network<sequential> net;

    // first input dicom image size
    // final output 2 classes

    // add layers
    /* Conv_layer
    * @param in_width         [in] input image width
    * @param in_height        [in] input image height
    * @param window_size      [in] window(kernel) size of convolution
    * @param in_channels      [in] input image channels (grayscale=1, rgb=3)
    * @param out_channels     [in] output image channels
    *   output-size = (in-width - window_width + 1) *
    * (in-height - window_height + 1) * out_channels
    */

    /* Average_pool layer
     * @param in_width     [in] width of input image
     * @param in_height    [in] height of input image
     * @param in_channels  [in] the number of input image channels(depth)
     * @param pool_size    [in] factor by which to downscale
     */
//    net << conv(32, 32, 5, 1, 6) << activation::tanh()  // in:32x32x1, 5x5conv, 6fmaps
//        << ave_pool(28, 28, 6, 2) << activation::tanh() // in:28x28x6, 2x2pooling
//        << fc(14 * 14 * 6, 120) << activation::tanh()   // in:14x14x6, out:120
//        << fc(120, 10);                     // in:120,     out:10

vec_t flatImage(vector<vector<double>> image) {

    vec_t flattedImage;
    for (auto row : image) {
        for (auto val : row) {
            flattedImage.push_back(val);
        }
    }

    return flattedImage;
}


int main(int argc, char *argv[])
{
    DicomReader dicomObj("/home/will/Projects/dicom-classifier/data/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");
    int width = dicomObj.getHeight();
    int height = dicomObj.getWidth();

    vector<vector<double>> imageData = dicomObj.getDoubleImageMatrix(12);
    vec_t flatImg = flatImage(imageData);

    vector<vec_t> input_data;
    vector<label_t> labels {0};

    input_data.push_back(flatImg);

    network<sequential> net;
    adam opt;

//    // add layers
    cout << "Adding layers... "<<endl;

    net << conv(width, height, 5, 1, 1) << activation::tanh()  // in: widthxheightx1, 5x5conv,
        << ave_pool(width-5+1, height-5+1, 1, 2) << activation::tanh() // in: widthxheightx1, 2x2pooling
        << fc(1662 * 2040 * 1, 10) << activation::tanh()   // in: 24x30, out:10
        << fc(10, 1);; //in: 10 out: 3 classes

    size_t batch_size = 1;
    size_t epochs = 1 ;

    cout << "Begin training"<<endl;
    net.train<mse>(opt, input_data, labels, batch_size, epochs);

//    net.save("weights", content_type::weights, file_format::json);

    for (int i = 0; i < net.depth(); i++) {
        cout << "#layer:" << i << "\n";
        cout << "layer type:" << net[i]->layer_type() << "\n";
        cout << "input:" << net[i]->in_data_size() << "(" << net[i]->in_data_shape() << ")\n";
        cout << "output:" << net[i]->out_data_size() << "(" << net[i]->out_data_shape() << ")\n";
    }


//    string json = net.to_json();
//    cout << json << endl;


    cout << "finished" << endl;
    return 0;
}

