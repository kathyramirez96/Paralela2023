/*
    Copyright (c) 2013, Taiga Nomi and the respective contributors
    All rights reserved.

    Use of this source code is governed by a BSD-style license that can be found
    in the LICENSE file.
*/
#pragma once

#include <string>
#include <utility>
#include <vector>

#include <neural_net/cereal/access.hpp>  // For LoadAndConstruct
#include "neural_net/nn.h"

namespace detail {

/**
 * size of layer, model, data etc.
 * change to smaller type if memory footprint is severe
 **/
typedef std::uint32_t serial_size_t;

typedef neural_net::index3d<serial_size_t> shape3d_serial;

template <class T>
static inline cereal::NameValuePair<T> make_nvp(const char *name, T &&value) {
  return cereal::make_nvp(name, value);
}

template <typename T>
struct is_binary_input_archive {
  static const bool value = false;
};
template <typename T>
struct is_binary_output_archive {
  static const bool value = false;
};
template <>
struct is_binary_input_archive<cereal::BinaryInputArchive> {
  static const bool value = true;
};
template <>
struct is_binary_input_archive<cereal::PortableBinaryInputArchive> {
  static const bool value = true;
};
template <>
struct is_binary_output_archive<cereal::BinaryOutputArchive> {
  static const bool value = true;
};
template <>
struct is_binary_output_archive<cereal::PortableBinaryOutputArchive> {
  static const bool value = true;
};

template <class Archive, typename dummy = Archive>
struct ArchiveWrapper {
  explicit ArchiveWrapper(Archive &ar) : ar(ar) {}
  template <typename T>
  void operator()(T &arg) {
    ar(arg);
  }
  Archive &ar;
};

template <typename Archive>
struct ArchiveWrapper<
  Archive,
  typename std::enable_if<is_binary_input_archive<Archive>::value,
                          Archive>::type> {
  explicit ArchiveWrapper(Archive &ar) : ar(ar) {}
  template <typename T>
  void operator()(T &arg) {
    ar(arg);
  }
  void operator()(cereal::NameValuePair<size_t &> &arg) {
    cereal::NameValuePair<serial_size_t> arg2(arg.name, 0);
    ar(arg2);
    arg.value = arg2.value;
  }
  Archive &ar;
};

template <typename Archive>
struct ArchiveWrapper<
  Archive,
  typename std::enable_if<is_binary_output_archive<Archive>::value,
                          Archive>::type> {
  explicit ArchiveWrapper(Archive &ar) : ar(ar) {}
  template <typename T>
  void operator()(T &arg) {
    ar(arg);
  }
  void operator()(cereal::NameValuePair<size_t &> &arg) {
    cereal::NameValuePair<serial_size_t> arg2(arg.name, 0);
    arg2.value = static_cast<serial_size_t>(arg.value);
    ar(arg2);
  }
  Archive &ar;
};

template <class Archive, typename T>
void arc(Archive &ar, T &&arg) {
  ArchiveWrapper<Archive> wa(ar);
  wa(arg);
}

template <class Archive>
inline void arc(Archive &ar) {}

template <class Archive, class Type, class Type2>
inline void arc(Archive &ar, Type &&arg, Type2 &&arg2) {
  arc(ar, std::forward<Type>(arg));
  arc(ar, std::forward<Type2>(arg2));
}

template <class Archive, class Type, class... Types>
inline void arc(Archive &ar, Type &&arg, Types &&... args) {
  arc(ar, std::forward<Type>(arg));
  arc(ar, std::forward<Types>(args)...);
}

}  // namespace detail

namespace cereal {

template <>
struct LoadAndConstruct<neural_net::elementwise_add_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::elementwise_add_layer> &construct) {
    size_t num_args, dim;

    ::detail::arc(ar, ::detail::make_nvp("num_args", num_args),
                  ::detail::make_nvp("dim", dim));
    construct(num_args, dim);
  }
};

template <>
struct LoadAndConstruct<neural_net::average_pooling_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::average_pooling_layer> &construct) {
    neural_net::shape3d in;
    size_t stride_x, stride_y, pool_size_x, pool_size_y;
    bool ceil_mode;
    neural_net::padding pad_type;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("pool_size_x", pool_size_x),
                  ::detail::make_nvp("pool_size_y", pool_size_y),
                  ::detail::make_nvp("stride_x", stride_x),
                  ::detail::make_nvp("stride_y", stride_y),
                  ::detail::make_nvp("ceil_mode", ceil_mode),
                  ::detail::make_nvp("pad_type", pad_type));
    construct(in.width_, in.height_, in.depth_, pool_size_x, pool_size_y,
              stride_x, stride_y, ceil_mode, pad_type);
  }
};

template <>
struct LoadAndConstruct<neural_net::average_unpooling_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::average_unpooling_layer> &construct) {
    neural_net::shape3d in;
    size_t pool_size, stride;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("pool_size", pool_size),
                  ::detail::make_nvp("stride", stride));
    construct(in.width_, in.height_, in.depth_, pool_size, stride);
  }
};

template <>
struct LoadAndConstruct<neural_net::batch_normalization_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::batch_normalization_layer> &construct) {
    size_t in_spatial_size, in_channels;
    neural_net::float_t eps, momentum;
    neural_net::net_phase phase;
    neural_net::vec_t mean, variance;

    ::detail::arc(ar, ::detail::make_nvp("in_spatial_size", in_spatial_size),
                  ::detail::make_nvp("in_channels", in_channels),
                  ::detail::make_nvp("epsilon", eps),
                  ::detail::make_nvp("momentum", momentum),
                  ::detail::make_nvp("phase", phase),
                  ::detail::make_nvp("mean", mean),
                  ::detail::make_nvp("variance", variance));
    construct(in_spatial_size, in_channels, eps, momentum, phase);
    construct->set_mean(mean);
    construct->set_variance(variance);
  }
};

template <>
struct LoadAndConstruct<neural_net::concat_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::concat_layer> &construct) {
    std::vector<neural_net::shape3d> in_shapes;
    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shapes));
    construct(in_shapes);
  }
};

template <>
struct LoadAndConstruct<neural_net::convolutional_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::convolutional_layer> &construct) {
    size_t w_width, w_height, out_ch, w_stride, h_stride, w_dilation,
      h_dilation;
    bool has_bias;
    neural_net::shape3d in;
    neural_net::padding pad_type;
    neural_net::core::connection_table tbl;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("window_width", w_width),
                  ::detail::make_nvp("window_height", w_height),
                  ::detail::make_nvp("out_channels", out_ch),
                  ::detail::make_nvp("connection_table", tbl),
                  ::detail::make_nvp("pad_type", pad_type),
                  ::detail::make_nvp("has_bias", has_bias),
                  ::detail::make_nvp("w_stride", w_stride),
                  ::detail::make_nvp("h_stride", h_stride),
                  ::detail::make_nvp("w_dilation", w_dilation),
                  ::detail::make_nvp("h_dilation", h_dilation));

    construct(in.width_, in.height_, w_width, w_height, in.depth_, out_ch, tbl,
              pad_type, has_bias, w_stride, h_stride, w_dilation, h_dilation);
  }
};

template <>
struct LoadAndConstruct<neural_net::deconvolutional_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::deconvolutional_layer> &construct) {
    size_t w_width, w_height, out_ch, w_stride, h_stride;
    bool has_bias;
    neural_net::shape3d in;
    neural_net::padding pad_type;
    neural_net::core::connection_table tbl;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("window_width", w_width),
                  ::detail::make_nvp("window_height", w_height),
                  ::detail::make_nvp("out_channels", out_ch),
                  ::detail::make_nvp("connection_table", tbl),
                  ::detail::make_nvp("pad_type", pad_type),
                  ::detail::make_nvp("has_bias", has_bias),
                  ::detail::make_nvp("w_stride", w_stride),
                  ::detail::make_nvp("h_stride", h_stride));

    construct(in.width_, in.height_, w_width, w_height, in.depth_, out_ch, tbl,
              pad_type, has_bias, w_stride, h_stride);
  }
};

template <>
struct LoadAndConstruct<neural_net::dropout_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::dropout_layer> &construct) {
    neural_net::net_phase phase;
    neural_net::float_t dropout_rate;
    size_t in_size;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_size),
                  ::detail::make_nvp("dropout_rate", dropout_rate),
                  ::detail::make_nvp("phase", phase));
    construct(in_size, dropout_rate, phase);
  }
};

template <>
struct LoadAndConstruct<neural_net::fully_connected_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::fully_connected_layer> &construct) {
    size_t in_dim, out_dim;
    bool has_bias;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_dim),
                  ::detail::make_nvp("out_size", out_dim),
                  ::detail::make_nvp("has_bias", has_bias));
    construct(in_dim, out_dim, has_bias);
  }
};

template <>
struct LoadAndConstruct<neural_net::global_average_pooling_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::global_average_pooling_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_shape", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::input_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::input_layer> &construct) {
    neural_net::shape3d shape;

    ::detail::arc(ar, ::detail::make_nvp("shape", shape));
    construct(shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::l2_normalization_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::l2_normalization_layer> &construct) {
    size_t in_spatial_size, in_channels;
    neural_net::float_t eps, scale;

    ::detail::arc(ar, ::detail::make_nvp("in_spatial_size", in_spatial_size),
                  ::detail::make_nvp("in_channels", in_channels),
                  ::detail::make_nvp("epsilon", eps),
                  ::detail::make_nvp("scale", scale));
    construct(in_spatial_size, in_channels, eps, scale);
  }
};

template <>
struct LoadAndConstruct<neural_net::linear_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::linear_layer> &construct) {
    size_t dim;
    neural_net::float_t scale, bias;

    ::detail::arc(ar, ::detail::make_nvp("in_size", dim),
                  ::detail::make_nvp("scale", scale),
                  ::detail::make_nvp("bias", bias));

    construct(dim, scale, bias);
  }
};

template <>
struct LoadAndConstruct<neural_net::lrn_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::lrn_layer> &construct) {
    neural_net::shape3d in_shape;
    size_t size;
    neural_net::float_t alpha, beta;
    neural_net::norm_region region;

    ::detail::arc(
      ar, ::detail::make_nvp("in_shape", in_shape),
      ::detail::make_nvp("size", size), ::detail::make_nvp("alpha", alpha),
      ::detail::make_nvp("beta", beta), ::detail::make_nvp("region", region));
    construct(in_shape, size, alpha, beta, region);
  }
};

template <>
struct LoadAndConstruct<neural_net::max_pooling_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::max_pooling_layer> &construct) {
    neural_net::shape3d in;
    size_t stride_x, stride_y, pool_size_x, pool_size_y;
    bool ceil_mode;
    neural_net::padding pad_type;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("pool_size_x", pool_size_x),
                  ::detail::make_nvp("pool_size_y", pool_size_y),
                  ::detail::make_nvp("stride_x", stride_x),
                  ::detail::make_nvp("stride_y", stride_y),
                  ::detail::make_nvp("ceil_mode", ceil_mode),
                  ::detail::make_nvp("pad_type", pad_type));
    construct(in.width_, in.height_, in.depth_, pool_size_x, pool_size_y,
              stride_x, stride_y, ceil_mode, pad_type);
  }
};

template <>
struct LoadAndConstruct<neural_net::max_unpooling_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::max_unpooling_layer> &construct) {
    neural_net::shape3d in;
    size_t stride, unpool_size;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("unpool_size", unpool_size),
                  ::detail::make_nvp("stride", stride));
    construct(in, unpool_size, stride);
  }
};

template <>
struct LoadAndConstruct<neural_net::zero_pad_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::zero_pad_layer> &construct) {
    neural_net::shape3d in_shape;
    size_t w_pad_size, h_pad_size;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("w_pad_size", w_pad_size),
                  ::detail::make_nvp("h_pad_size", h_pad_size));
    construct(in_shape, w_pad_size, h_pad_size);
  }
};

template <>
struct LoadAndConstruct<neural_net::power_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::power_layer> &construct) {
    neural_net::shape3d in_shape;
    neural_net::float_t factor;
    neural_net::float_t scale(1.0f);

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("factor", factor),
                  ::detail::make_nvp("scale", scale));
    construct(in_shape, factor, scale);
  }
};

template <>
struct LoadAndConstruct<neural_net::quantized_convolutional_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::quantized_convolutional_layer> &construct) {
    size_t w_width, w_height, out_ch, w_stride, h_stride;
    bool has_bias;
    neural_net::shape3d in;
    neural_net::padding pad_type;
    neural_net::core::connection_table tbl;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("window_width", w_width),
                  ::detail::make_nvp("window_height", w_height),
                  ::detail::make_nvp("out_channels", out_ch),
                  ::detail::make_nvp("connection_table", tbl),
                  ::detail::make_nvp("pad_type", pad_type),
                  ::detail::make_nvp("has_bias", has_bias),
                  ::detail::make_nvp("w_stride", w_stride),
                  ::detail::make_nvp("h_stride", h_stride));

    construct(in.width_, in.height_, w_width, w_height, in.depth_, out_ch, tbl,
              pad_type, has_bias, w_stride, h_stride);
  }
};

template <>
struct LoadAndConstruct<neural_net::quantized_deconvolutional_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::quantized_deconvolutional_layer> &construct) {
    size_t w_width, w_height, out_ch, w_stride, h_stride;
    bool has_bias;
    neural_net::shape3d in;
    neural_net::padding pad_type;
    neural_net::core::connection_table tbl;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in),
                  ::detail::make_nvp("window_width", w_width),
                  ::detail::make_nvp("window_height", w_height),
                  ::detail::make_nvp("out_channels", out_ch),
                  ::detail::make_nvp("connection_table", tbl),
                  ::detail::make_nvp("pad_type", pad_type),
                  ::detail::make_nvp("has_bias", has_bias),
                  ::detail::make_nvp("w_stride", w_stride),
                  ::detail::make_nvp("h_stride", h_stride));

    construct(in.width_, in.height_, w_width, w_height, in.depth_, out_ch, tbl,
              pad_type, has_bias, w_stride, h_stride);
  }
};

template <>
struct LoadAndConstruct<neural_net::quantized_fully_connected_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar,
    cereal::construct<neural_net::quantized_fully_connected_layer> &construct) {
    size_t in_dim, out_dim;
    bool has_bias;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_dim),
                  ::detail::make_nvp("out_size", out_dim),
                  ::detail::make_nvp("has_bias", has_bias));
    construct(in_dim, out_dim, has_bias);
  }
};

template <>
struct LoadAndConstruct<neural_net::recurrent_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::recurrent_layer> &construct) {
    size_t seq_len;
    ::detail::arc(ar, ::detail::make_nvp("seq_len", seq_len));
    auto cell_p = neural_net::layer::load_layer(ar);

    construct(std::static_pointer_cast<neural_net::cell>(cell_p), seq_len);
  }
};

template <>
struct LoadAndConstruct<neural_net::gru_cell> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::gru_cell> &construct) {
    size_t in_dim, out_dim;
    bool has_bias;
    ::detail::arc(ar, ::detail::make_nvp("in_size", in_dim),
                  ::detail::make_nvp("out_size", out_dim),
                  ::detail::make_nvp("has_bias", has_bias));
    neural_net::gru_cell_parameters params;
    params.has_bias = has_bias;
    construct(in_dim, out_dim, params);
  }
};

template <>
struct LoadAndConstruct<neural_net::lstm_cell> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::lstm_cell> &construct) {
    size_t in_dim, out_dim;
    bool has_bias;
    ::detail::arc(ar, ::detail::make_nvp("in_size", in_dim),
                  make_nvp("out_size", out_dim),
                  ::detail::make_nvp("has_bias", has_bias));
    neural_net::lstm_cell_parameters params;
    params.has_bias = has_bias;
    construct(in_dim, out_dim, params);
  }
};

template <>
struct LoadAndConstruct<neural_net::rnn_cell> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::rnn_cell> &construct) {
    size_t in_dim, out_dim;
    bool has_bias;
    ::detail::arc(ar, ::detail::make_nvp("in_size", in_dim),
                  ::detail::make_nvp("out_size", out_dim),
                  ::detail::make_nvp("has_bias", has_bias));
    neural_net::rnn_cell_parameters params;
    params.has_bias = has_bias;
    construct(in_dim, out_dim, params);
  }
};

template <>
struct LoadAndConstruct<neural_net::slice_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::slice_layer> &construct) {
    neural_net::shape3d in_shape;
    neural_net::slice_type slice_type;
    size_t num_outputs;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("slice_type", slice_type),
                  ::detail::make_nvp("num_outputs", num_outputs));
    construct(in_shape, slice_type, num_outputs);
  }
};

template <>
struct LoadAndConstruct<neural_net::sigmoid_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::sigmoid_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::asinh_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::asinh_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::tanh_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::tanh_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::relu_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::relu_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::softmax_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::softmax_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::leaky_relu_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::leaky_relu_layer> &construct) {
    neural_net::shape3d in_shape;
    neural_net::float_t epsilon;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("epsilon", epsilon));

    construct(in_shape, epsilon);
  }
};

template <>
struct LoadAndConstruct<neural_net::selu_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::selu_layer> &construct) {
    neural_net::shape3d in_shape;
    neural_net::float_t lambda;
    neural_net::float_t alpha;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("lambda", lambda),
                  ::detail::make_nvp("alpha", alpha));
    construct(in_shape, lambda, alpha);
  }
};

template <>
struct LoadAndConstruct<neural_net::elu_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::elu_layer> &construct) {
    neural_net::shape3d in_shape;
    neural_net::float_t alpha;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("alpha", alpha));
    construct(in_shape, alpha);
  }
};

template <>
struct LoadAndConstruct<neural_net::tanh_p1m2_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::tanh_p1m2_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

template <>
struct LoadAndConstruct<neural_net::softplus_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::softplus_layer> &construct) {
    neural_net::shape3d in_shape;
    neural_net::float_t beta;
    neural_net::float_t threshold;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape),
                  ::detail::make_nvp("beta", beta),
                  ::detail::make_nvp("threshold", threshold));
    construct(in_shape, beta, threshold);
  }
};

template <>
struct LoadAndConstruct<neural_net::softsign_layer> {
  template <class Archive>
  static void load_and_construct(
    Archive &ar, cereal::construct<neural_net::softsign_layer> &construct) {
    neural_net::shape3d in_shape;

    ::detail::arc(ar, ::detail::make_nvp("in_size", in_shape));
    construct(in_shape);
  }
};

}  // namespace cereal

namespace neural_net {

struct serialization_buddy {
#ifndef CNN_NO_SERIALIZATION

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::layer &layer) {
    auto all_weights = layer.weights();
    for (auto weight : all_weights) {
      ar(*weight);
    }
    layer.initialized_ = true;
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::elementwise_add_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("num_args", layer.num_args_),
                  ::detail::make_nvp("dim", layer.dim_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::average_pooling_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_),
                  ::detail::make_nvp("pool_size_x", layer.pool_size_x_),
                  ::detail::make_nvp("pool_size_y", layer.pool_size_y_),
                  ::detail::make_nvp("stride_x", layer.stride_x_),
                  ::detail::make_nvp("stride_y", layer.stride_y_),
                  ::detail::make_nvp("ceil_mode", layer.ceil_mode_),
                  ::detail::make_nvp("pad_type", layer.pad_type_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::average_unpooling_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_),
                  ::detail::make_nvp("pool_size", layer.w_.width_),
                  ::detail::make_nvp("stride", layer.stride_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::batch_normalization_layer &layer) {
    ::detail::arc(ar,
                  ::detail::make_nvp("in_spatial_size", layer.in_spatial_size_),
                  ::detail::make_nvp("in_channels", layer.in_channels_),
                  ::detail::make_nvp("epsilon", layer.eps_),
                  ::detail::make_nvp("momentum", layer.momentum_),
                  ::detail::make_nvp("phase", layer.phase_),
                  ::detail::make_nvp("mean", layer.mean_),
                  ::detail::make_nvp("variance", layer.variance_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::concat_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shapes_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::convolutional_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in),
                  ::detail::make_nvp("window_width", params_.weight.width_),
                  ::detail::make_nvp("window_height", params_.weight.height_),
                  ::detail::make_nvp("out_channels", params_.out.depth_),
                  ::detail::make_nvp("connection_table", params_.tbl),
                  ::detail::make_nvp("pad_type", params_.pad_type),
                  ::detail::make_nvp("has_bias", params_.has_bias),
                  ::detail::make_nvp("w_stride", params_.w_stride),
                  ::detail::make_nvp("h_stride", params_.w_stride),
                  ::detail::make_nvp("w_dilation", params_.w_dilation),
                  ::detail::make_nvp("h_dilation", params_.h_dilation));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::deconvolutional_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in),
                  ::detail::make_nvp("window_width", params_.weight.width_),
                  ::detail::make_nvp("window_height", params_.weight.height_),
                  ::detail::make_nvp("out_channels", params_.out.depth_),
                  ::detail::make_nvp("connection_table", params_.tbl),
                  ::detail::make_nvp("pad_type", params_.pad_type),
                  ::detail::make_nvp("has_bias", params_.has_bias),
                  ::detail::make_nvp("w_stride", params_.w_stride),
                  ::detail::make_nvp("h_stride", params_.h_stride));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::dropout_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_size_),
                  ::detail::make_nvp("dropout_rate", layer.dropout_rate_),
                  ::detail::make_nvp("phase", layer.phase_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::fully_connected_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in_size_),
                  ::detail::make_nvp("out_size", params_.out_size_),
                  ::detail::make_nvp("has_bias", params_.has_bias_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::global_average_pooling_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_shape", params_.in));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::input_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("shape", layer.shape_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::l2_normalization_layer &layer) {
    ::detail::arc(ar,
                  ::detail::make_nvp("in_spatial_size", layer.in_spatial_size_),
                  ::detail::make_nvp("in_channels", layer.in_channels_),
                  ::detail::make_nvp("epsilon", layer.eps_),
                  ::detail::make_nvp("scale", layer.scale_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::linear_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.dim_),
                  ::detail::make_nvp("scale", layer.scale_),
                  ::detail::make_nvp("bias", layer.bias_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::lrn_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_shape", layer.in_shape_),
                  ::detail::make_nvp("size", layer.size_),
                  ::detail::make_nvp("alpha", layer.alpha_),
                  ::detail::make_nvp("beta", layer.beta_),
                  ::detail::make_nvp("region", layer.region_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::max_pooling_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in),
                  ::detail::make_nvp("pool_size_x", params_.pool_size_x),
                  ::detail::make_nvp("pool_size_y", params_.pool_size_y),
                  ::detail::make_nvp("stride_x", params_.stride_x),
                  ::detail::make_nvp("stride_y", params_.stride_y),
                  ::detail::make_nvp("ceil_mode", params_.ceil_mode),
                  ::detail::make_nvp("pad_type", params_.pad_type));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::max_unpooling_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_),
                  ::detail::make_nvp("unpool_size", layer.unpool_size_),
                  ::detail::make_nvp("stride", layer.stride_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::zero_pad_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape_),
                  ::detail::make_nvp("w_pad_size", layer.w_pad_size_),
                  ::detail::make_nvp("h_pad_size", layer.h_pad_size_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::power_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape_),
                  ::detail::make_nvp("factor", layer.factor_),
                  ::detail::make_nvp("scale", layer.scale_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar,
                               neural_net::quantized_convolutional_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in),
                  ::detail::make_nvp("window_width", params_.weight.width_),
                  ::detail::make_nvp("window_height", params_.weight.height_),
                  ::detail::make_nvp("out_channels", params_.out.depth_),
                  ::detail::make_nvp("connection_table", params_.tbl),
                  ::detail::make_nvp("pad_type", params_.pad_type),
                  ::detail::make_nvp("has_bias", params_.has_bias),
                  ::detail::make_nvp("w_stride", params_.w_stride),
                  ::detail::make_nvp("h_stride", params_.h_stride));
  }

  template <class Archive>
  static inline void serialize(
    Archive &ar, neural_net::quantized_deconvolutional_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in),
                  ::detail::make_nvp("window_width", params_.weight.width_),
                  ::detail::make_nvp("window_height", params_.weight.height_),
                  ::detail::make_nvp("out_channels", params_.out.depth_),
                  ::detail::make_nvp("connection_table", params_.tbl),
                  ::detail::make_nvp("pad_type", params_.pad_type),
                  ::detail::make_nvp("has_bias", params_.has_bias),
                  ::detail::make_nvp("w_stride", params_.w_stride),
                  ::detail::make_nvp("h_stride", params_.h_stride));
  }

  template <class Archive>
  static inline void serialize(
    Archive &ar, neural_net::quantized_fully_connected_layer &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in_size_),
                  ::detail::make_nvp("out_size", params_.out_size_),
                  ::detail::make_nvp("has_bias", params_.has_bias_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::recurrent_layer &layer) {
    size_t seq_len = layer.seq_len_;
    ::detail::arc(ar, ::detail::make_nvp("seq_len", seq_len));
    neural_net::layer::save_layer(ar, *layer.cell_);
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::gru_cell &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in_size_),
                  ::detail::make_nvp("out_size", params_.out_size_),
                  ::detail::make_nvp("has_bias", params_.has_bias_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::lstm_cell &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in_size_),
                  ::detail::make_nvp("out_size", params_.out_size_),
                  ::detail::make_nvp("has_bias", params_.has_bias_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::rnn_cell &layer) {
    auto &params_ = layer.params_;
    ::detail::arc(ar, ::detail::make_nvp("in_size", params_.in_size_),
                  ::detail::make_nvp("out_size", params_.out_size_),
                  ::detail::make_nvp("has_bias", params_.has_bias_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::slice_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape_),
                  ::detail::make_nvp("slice_type", layer.slice_type_),
                  ::detail::make_nvp("num_outputs", layer.num_outputs_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::sigmoid_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::asinh_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::tanh_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::relu_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::softmax_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::leaky_relu_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]),
                  ::detail::make_nvp("epsilon", layer.epsilon_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::elu_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]),
                  ::detail::make_nvp("alpha", layer.alpha_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::selu_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]),
                  ::detail::make_nvp("lambda", layer.lambda_),
                  ::detail::make_nvp("alpha", layer.alpha_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::tanh_p1m2_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::softplus_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]),
                  ::detail::make_nvp("beta", layer.beta_),
                  ::detail::make_nvp("threshold", layer.threshold_));
  }

  template <class Archive>
  static inline void serialize(Archive &ar, neural_net::softsign_layer &layer) {
    ::detail::arc(ar, ::detail::make_nvp("in_size", layer.in_shape()[0]));
  }

#endif  // #ifndef CNN_NO_SERIALIZATION
};      // struct serialization_buddy

template <class Archive, typename T>
typename std::enable_if<std::is_base_of<neural_net::layer, T>::value>::type
serialize(Archive &ar, T &layer) {
  auto &inconstant_layer =
    const_cast<typename std::remove_const<T>::type &>(layer);
  inconstant_layer.serialize_prolog(ar);
  serialization_buddy::serialize(ar, inconstant_layer);
}

template <class Archive, typename T>
void serialize(Archive &ar, neural_net::index3d<T> &idx) {
  ::detail::arc(ar, ::detail::make_nvp("width", idx.width_),
                ::detail::make_nvp("height", idx.height_),
                ::detail::make_nvp("depth", idx.depth_));
}

namespace core {

template <class Archive>
void serialize(Archive &ar, neural_net::core::connection_table &tbl) {
  ::detail::arc(ar, ::detail::make_nvp("rows", tbl.rows_),
                ::detail::make_nvp("cols", tbl.cols_));
  if (tbl.is_empty()) {
    std::string all("all");
    ::detail::arc(ar, ::detail::make_nvp("connection", all));
  } else {
    ::detail::arc(ar, ::detail::make_nvp("connection", tbl.connected_));
  }
}

}  // namespace core

}  // namespace neural_net
