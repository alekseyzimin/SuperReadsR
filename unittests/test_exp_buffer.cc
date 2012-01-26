#include <gtest/gtest.h>
#include <exp_buffer.hpp>
#include <algorithm>

template<typename T>
class ExpandingBufferInit : public ::testing::Test { };
template<typename T>
class ExpandingBufferDefault : public ::testing::Test { };

typedef ExpandingBuffer<int, reallocator<int> > int_realloc_buf;
typedef ExpandingBuffer<int, reallocator_init<int, -1> > init_realloc_buf;
typedef ExpandingBuffer<int, remaper<int> > int_remap_buf;
typedef ExpandingBuffer<int, remaper_init<int, -1> > init_remap_buf;

typedef ::testing::Types<int_realloc_buf, int_remap_buf> ExpBufferDefaultTypes;
typedef ::testing::Types<init_realloc_buf, init_remap_buf> ExpBufferInitTypes;

TYPED_TEST_CASE(ExpandingBufferDefault, ExpBufferDefaultTypes);
TYPED_TEST_CASE(ExpandingBufferInit, ExpBufferInitTypes);

TYPED_TEST(ExpandingBufferDefault, Initialization) {
  TypeParam b;

  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.size());

  b[5] = 5;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.size());
  // This holds only for remaper
  // for(int i = 0; i < 5; ++i)
  //   EXPECT_EQ(0, b[i]);
  EXPECT_EQ(5, b[5]);
  b[3] = 3;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.size());
  EXPECT_EQ(3, b[3]);

  b[6] = 6;
  EXPECT_EQ((size_t)12, b.capacity());
  EXPECT_EQ((size_t)7, b.size());

  b[5000] = 5000;
  EXPECT_EQ((size_t)5001, b.size());
  EXPECT_EQ(5000, b.back());
  // Again, only for remaper
  // for(typename TypeParam::iterator it = b.begin(); it != b.end(); ++it)
  //   EXPECT_TRUE(*it == 0 || *it == (it - b.begin()));

  b.push_back(5001);
  EXPECT_EQ((size_t)5002, b.size());
  EXPECT_EQ(5001, b.back());
}

TYPED_TEST(ExpandingBufferDefault, Swap) {
  TypeParam b(10);

  for(size_t i = 0; i < b.capacity(); ++i)
    b[i] = 2 * i;
  
  EXPECT_EQ((size_t)10, b.size());
  
  TypeParam bs;
  b.swap(bs);
  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)10, bs.capacity());
  for(size_t i = 0; i < bs.capacity(); ++i)
    EXPECT_EQ((int)(2 * i), bs[i]);

  std::swap(b, bs);
  
  EXPECT_EQ((size_t)10, b.capacity());
  EXPECT_EQ((size_t)0, bs.capacity());
  for(size_t i = 0; i < b.capacity(); ++i)
    EXPECT_EQ((int)(2 * i), b[i]);

  ExpBuffer<int> b1, b2(5);
  std::swap(b1, b2);
  EXPECT_EQ((size_t)5, b1.capacity());
  EXPECT_EQ((size_t)0, b2.capacity());
}


TYPED_TEST(ExpandingBufferInit, Initialization) {
  TypeParam b;
  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.size());
  b[10] = 5;
  EXPECT_EQ(b.size() - 1, (size_t)std::count(b.begin(), b.end(), -1));
}

