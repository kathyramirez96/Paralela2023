QT -= gui

TEMPLATE = lib
CONFIG += staticlib

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

QMAKE_CXXFLAGS += -pthread

LIBS += -lpthread



SOURCES += \
    dicom_read/DicomReader.cpp \
    dicom_read/dicomutils.cpp \
    knn/Preprocessing.cpp \
    knn/ReadDataset.cpp \
    knn/dataset.cpp \
    knn/knn.cpp \
    knn/knnresults.cpp \
    knn/matrix.cpp \
    random_forest/ArgumentHandler.cpp \
    random_forest/Data.cpp \
    random_forest/Forest.cpp \
    random_forest/ForestClassification.cpp \
    random_forest/ForestProbability.cpp \
    random_forest/ForestRegression.cpp \
    random_forest/ForestSurvival.cpp \
    random_forest/Tree.cpp \
    random_forest/TreeClassification.cpp \
    random_forest/TreeProbability.cpp \
    random_forest/TreeRegression.cpp \
    random_forest/TreeSurvival.cpp \
    random_forest/rf.cpp \
    random_forest/utility.cpp \
    svm/svm.cpp \
    fcm/fcm.cpp \
    svm/svm_util.cpp


HEADERS += \
    dicom_read/DicomReader.h \
    dicom_read/dicomutils.h \
    knn/Preprocessing.h \
    knn/ReadDataset.h \
    knn/dataset.h \
    knn/debug.h \
    knn/knn.h \
    knn/knnUtils.h \
    knn/matrix.h \
    random_forest/rf.h \
    svm/svm.h \
    fcm/fcm.h \
    kmeans/kmeans.h \
    svm/svm_util.h



INCLUDEPATH += random_forest

# Default rules for deployment.
unix {
    target.path = $$[QT_INSTALL_PLUGINS]/generic
}
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    tiny_dnn/io/caffe/CPPLINT.cfg \
    tiny_dnn/io/caffe/caffe.proto \
    tiny_dnn/third_party/cl_kernels/activation.cl \
    tiny_dnn/third_party/cl_kernels/auxiliary.cl \
    tiny_dnn/third_party/cl_kernels/batch_reindex.cl \
    tiny_dnn/third_party/cl_kernels/benchmark.cl \
    tiny_dnn/third_party/cl_kernels/bias.cl \
    tiny_dnn/third_party/cl_kernels/bnll.cl \
    tiny_dnn/third_party/cl_kernels/channel.cl \
    tiny_dnn/third_party/cl_kernels/concat.cl \
    tiny_dnn/third_party/cl_kernels/contrastive_loss.cl \
    tiny_dnn/third_party/cl_kernels/conv_layer_spatial.cl \
    tiny_dnn/third_party/cl_kernels/conv_spatial_helper.cl \
    tiny_dnn/third_party/cl_kernels/crop.cl \
    tiny_dnn/third_party/cl_kernels/dropout.cl \
    tiny_dnn/third_party/cl_kernels/eltwise.cl \
    tiny_dnn/third_party/cl_kernels/elu.cl \
    tiny_dnn/third_party/cl_kernels/embed.cl \
    tiny_dnn/third_party/cl_kernels/fft.cl \
    tiny_dnn/third_party/cl_kernels/fillbuffer.cl \
    tiny_dnn/third_party/cl_kernels/im2col.cl \
    tiny_dnn/third_party/cl_kernels/im2col_nd.cl \
    tiny_dnn/third_party/cl_kernels/lrn.cl \
    tiny_dnn/third_party/cl_kernels/lstm_unit.cl \
    tiny_dnn/third_party/cl_kernels/math.cl \
    tiny_dnn/third_party/cl_kernels/mergecrop.cl \
    tiny_dnn/third_party/cl_kernels/pooling.cl \
    tiny_dnn/third_party/cl_kernels/pooling_nd.cl \
    tiny_dnn/third_party/cl_kernels/pooling_sk.cl \
    tiny_dnn/third_party/cl_kernels/slice.cl \
    tiny_dnn/third_party/cl_kernels/softmax_loss.cl \
    tiny_dnn/third_party/cl_kernels/solvers.cl \
    tiny_dnn/third_party/cl_kernels/tile.cl \
    tiny_dnn/third_party/cpplint.py \
    tiny_dnn/third_party/gemmlowp/meta/README \
    tiny_dnn/third_party/gemmlowp/meta/generators/cc_emitter.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/gemm_MxNxK.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/gemv_1xMxK.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/meta.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/meta_arm_32.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/meta_arm_64.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/meta_arm_common.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/mul_1x8_Mx8_neon.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/mul_Nx8_Mx8_neon.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/neon_emitter.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/neon_emitter_64.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/qnt_Nx8_neon.py \
    tiny_dnn/third_party/gemmlowp/meta/generators/zip_Nx8_neon.py
