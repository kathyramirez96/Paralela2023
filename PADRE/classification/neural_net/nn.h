/*
    Copyright (c) 2013, Taiga Nomi and the respective contributors
    All rights reserved.

    Use of this source code is governed by a BSD-style license that can be found
    in the LICENSE file.
*/
#pragma once

#include "neural_net/config.h"
#include "neural_net/network.h"
#include "neural_net/nodes.h"

#include "neural_net/core/framework/tensor.h"

#include "neural_net/core/framework/device.h"
#include "neural_net/core/framework/program_manager.h"

#include "neural_net/activations/asinh_layer.h"
#include "neural_net/activations/elu_layer.h"
#include "neural_net/activations/leaky_relu_layer.h"
#include "neural_net/activations/relu_layer.h"
#include "neural_net/activations/selu_layer.h"
#include "neural_net/activations/sigmoid_layer.h"
#include "neural_net/activations/softmax_layer.h"
#include "neural_net/activations/softplus_layer.h"
#include "neural_net/activations/softsign_layer.h"
#include "neural_net/activations/tanh_layer.h"
#include "neural_net/activations/tanh_p1m2_layer.h"
#include "neural_net/layers/arithmetic_layer.h"
#include "neural_net/layers/average_pooling_layer.h"
#include "neural_net/layers/average_unpooling_layer.h"
#include "neural_net/layers/batch_normalization_layer.h"
#include "neural_net/layers/cell.h"
#include "neural_net/layers/cells.h"
#include "neural_net/layers/concat_layer.h"
#include "neural_net/layers/convolutional_layer.h"
#include "neural_net/layers/deconvolutional_layer.h"
#include "neural_net/layers/dropout_layer.h"
#include "neural_net/layers/fully_connected_layer.h"
#include "neural_net/layers/global_average_pooling_layer.h"
#include "neural_net/layers/input_layer.h"
#include "neural_net/layers/l2_normalization_layer.h"
#include "neural_net/layers/lrn_layer.h"
#include "neural_net/layers/lrn_layer.h"
#include "neural_net/layers/max_pooling_layer.h"
#include "neural_net/layers/max_unpooling_layer.h"
#include "neural_net/layers/power_layer.h"
#include "neural_net/layers/quantized_convolutional_layer.h"
#include "neural_net/layers/quantized_deconvolutional_layer.h"
#include "neural_net/layers/recurrent_layer.h"
#include "neural_net/layers/slice_layer.h"
#include "neural_net/layers/zero_pad_layer.h"

#ifdef CNN_USE_GEMMLOWP
#include "neural_net/layers/quantized_fully_connected_layer.h"
#endif  // CNN_USE_GEMMLOWP

#include "neural_net/lossfunctions/loss_function.h"
#include "neural_net/optimizers/optimizer.h"

#include "neural_net/util/deform.h"
#include "neural_net/util/graph_visualizer.h"
#include "neural_net/util/product.h"
#include "neural_net/util/weight_init.h"
#include "neural_net/util/nms.h"

#include "neural_net/io/cifar10_parser.h"
#include "neural_net/io/display.h"
#include "neural_net/io/layer_factory.h"
#include "neural_net/io/mnist_parser.h"

#ifdef DNN_USE_IMAGE_API
#include "neural_net/util/image.h"
#endif  // DNN_USE_IMAGE_API

#ifndef CNN_NO_SERIALIZATION
#include "neural_net/util/deserialization_helper.h"
#include "neural_net/util/serialization_helper.h"
// to allow upcasting
CEREAL_REGISTER_TYPE(neural_net::elu_layer)
CEREAL_REGISTER_TYPE(neural_net::leaky_relu_layer)
CEREAL_REGISTER_TYPE(neural_net::relu_layer)
CEREAL_REGISTER_TYPE(neural_net::sigmoid_layer)
CEREAL_REGISTER_TYPE(neural_net::softmax_layer)
CEREAL_REGISTER_TYPE(neural_net::softplus_layer)
CEREAL_REGISTER_TYPE(neural_net::softsign_layer)
CEREAL_REGISTER_TYPE(neural_net::tanh_layer)
CEREAL_REGISTER_TYPE(neural_net::tanh_p1m2_layer)
#endif  // CNN_NO_SERIALIZATION

// shortcut version of layer names
namespace neural_net {
namespace layers {

using conv = neural_net::convolutional_layer;

using q_conv = neural_net::quantized_convolutional_layer;

using max_pool = neural_net::max_pooling_layer;

using ave_pool = neural_net::average_pooling_layer;

using fc = neural_net::fully_connected_layer;

using dense = neural_net::fully_connected_layer;

using zero_pad = neural_net::zero_pad_layer;

// using rnn_cell = neural_net::rnn_cell_layer;

#ifdef CNN_USE_GEMMLOWP
using q_fc = neural_net::quantized_fully_connected_layer;
#endif

using add = neural_net::elementwise_add_layer;

using dropout = neural_net::dropout_layer;

using input = neural_net::input_layer;

using linear = linear_layer;

using lrn = neural_net::lrn_layer;

using concat = neural_net::concat_layer;

using deconv = neural_net::deconvolutional_layer;

using max_unpool = neural_net::max_unpooling_layer;

using ave_unpool = neural_net::average_unpooling_layer;

}  // namespace layers

namespace activation {

using sigmoid = neural_net::sigmoid_layer;

using asinh = neural_net::asinh_layer;

using tanh = neural_net::tanh_layer;

using relu = neural_net::relu_layer;

using rectified_linear = neural_net::relu_layer;

using softmax = neural_net::softmax_layer;

using leaky_relu = neural_net::leaky_relu_layer;

using elu = neural_net::elu_layer;

using selu = neural_net::selu_layer;

using tanh_p1m2 = neural_net::tanh_p1m2_layer;

using softplus = neural_net::softplus_layer;

using softsign = neural_net::softsign_layer;

}  // namespace activation

#include "neural_net/models/alexnet.h"

using batch_norm = neural_net::batch_normalization_layer;

using l2_norm = neural_net::l2_normalization_layer;

using slice = neural_net::slice_layer;

using power = neural_net::power_layer;

}  // namespace neural_net

#ifdef CNN_USE_CAFFE_CONVERTER
// experimental / require google protobuf
#include "neural_net/io/caffe/layer_factory.h"
#endif
