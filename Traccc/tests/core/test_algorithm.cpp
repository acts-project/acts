/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/utils/algorithm.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

class double_int : public traccc::algorithm<int(const int &)> {
    public:
    virtual int operator()(const int &i) const override { return 2 * i; }
};

class double_int_affine : public traccc::algorithm<int(int &&)> {
    public:
    virtual int operator()(int &&i) const override { return 2 * i; }
};

class double_string_affine
    : public traccc::algorithm<std::string(std::string &&)> {
    public:
    virtual std::string operator()(std::string &&i) const override {
        return i + i;
    }
};

class double_string_regular
    : public traccc::algorithm<std::string(const std::string &)> {
    public:
    virtual std::string operator()(const std::string &i) const override {
        return i + i;
    }
};

class add : public traccc::algorithm<int(const int &, const int &)> {
    public:
    virtual int operator()(const int &i, const int &j) const override {
        return i + j;
    }
};

TEST(algorithm, basic) {
    double_int i;

    ASSERT_EQ(i(1), 2);
    ASSERT_EQ(i(-1), -2);
    ASSERT_EQ(i(5), 10);
}

TEST(algorithm, basic_affine) {
    double_int_affine i;

    ASSERT_EQ(i(1), 2);
    ASSERT_EQ(i(-1), -2);
    ASSERT_EQ(i(5), 10);
}

TEST(algorithm, string_affine) {
    double_string_affine i;

    std::string s = "test";

    ASSERT_EQ(i("hello"), "hellohello");
    ASSERT_EQ(i("bye"), "byebye");
    ASSERT_EQ(i(std::move(s)), "testtest");
}

TEST(algorithm, string_regular) {
    double_string_regular i;

    std::string s = "test";

    ASSERT_EQ(i("hello"), "hellohello");
    ASSERT_EQ(i("bye"), "byebye");
    ASSERT_EQ(i(s), "testtest");
}

TEST(algorithm, compose_double_int_affine_1) {
    double_int_affine i;

    std::function<int(int &&)> f = traccc::compose(
        std::function<int(int &&)>(i), std::function<int(int &&)>(i));

    ASSERT_EQ(f(1), 4);
    ASSERT_EQ(f(-1), -4);
    ASSERT_EQ(f(5), 20);
}

TEST(algorithm, compose_double_int_affine_2) {
    double_int_affine i;
    double_int_affine j;

    std::function<int(int &&)> f = traccc::compose(
        std::function<int(int &&)>(i), std::function<int(int &&)>(j));

    ASSERT_EQ(f(1), 4);
    ASSERT_EQ(f(-1), -4);
    ASSERT_EQ(f(5), 20);
}

TEST(algorithm, compose_double_int_regular_1) {
    double_int i;

    std::function<int(const int &)> f = traccc::compose(
        std::function<int(const int &)>(i), std::function<int(const int &)>(i));

    ASSERT_EQ(f(1), 4);
    ASSERT_EQ(f(-1), -4);
    ASSERT_EQ(f(5), 20);
}

TEST(algorithm, compose_double_int_regular_2) {
    double_int i;
    double_int j;

    std::function<int(const int &)> f = traccc::compose(
        std::function<int(const int &)>(i), std::function<int(const int &)>(j));

    ASSERT_EQ(f(1), 4);
    ASSERT_EQ(f(-1), -4);
    ASSERT_EQ(f(5), 20);
}

TEST(algorithm, compose_string_lvalue_1) {
    double_string_regular i;

    std::string s = "test";

    std::function<std::string(const std::string &)> f =
        traccc::compose(std::function<std::string(const std::string &)>(i),
                        std::function<std::string(const std::string &)>(i));

    ASSERT_EQ(f("hello"), "hellohellohellohello");
    ASSERT_EQ(f("bye"), "byebyebyebye");
    ASSERT_EQ(f(s), "testtesttesttest");
}

TEST(algorithm, compose_string_rvalue_1) {
    double_string_affine i;

    std::string s = "test";

    std::function<std::string(std::string &&)> f =
        traccc::compose(std::function<std::string(std::string &&)>(i),
                        std::function<std::string(std::string &&)>(i));

    ASSERT_EQ(f("hello"), "hellohellohellohello");
    ASSERT_EQ(f("bye"), "byebyebyebye");
    ASSERT_EQ(f(std::move(s)), "testtesttesttest");
}

TEST(algorithm, compose_double_int_many_times) {
    double_int i;

    std::function<int(const int &)> f = traccc::compose(
        std::function<int(const int &)>(i), std::function<int(const int &)>(i),
        std::function<int(const int &)>(i), std::function<int(const int &)>(i),
        std::function<int(const int &)>(i), std::function<int(const int &)>(i),
        std::function<int(const int &)>(i), std::function<int(const int &)>(i));

    ASSERT_EQ(f(1), 256);
    ASSERT_EQ(f(-1), -256);
    ASSERT_EQ(f(5), 1280);
}

TEST(algorithm, side_effect) {
    double_int i;
    int j = 0;

    std::function<int(const int &)> f = traccc::compose(
        std::function<int(const int &)>(i), std::function<int(const int &)>(i),
        traccc::side_effect(std::function<int(const int &)>(i),
                            std::function<void(const int &)>(
                                [&j](const int &q) { j = q + 5; })));

    ASSERT_EQ(f(1), 8);
    ASSERT_EQ(j, 9);
}

TEST(algorithm, partial_application_bind) {
    add i;
    std::function<int(const int &)> f = std::bind(i, std::placeholders::_1, 5);

    ASSERT_EQ(f(1), 6);
    ASSERT_EQ(f(-1), 4);
    ASSERT_EQ(f(5), 10);
}

TEST(algorithm, partial_application_lambda) {
    add i;
    std::function<int(const int &)> f = [i](const int &j) { return i(j, 5); };

    ASSERT_EQ(f(1), 6);
    ASSERT_EQ(f(-1), 4);
    ASSERT_EQ(f(5), 10);
}

TEST(algorithm, partial_application_compose) {
    add i;
    std::function<int(const int &)> j = std::bind(i, std::placeholders::_1, 5);

    std::function<int(const int &)> f = traccc::compose(j, j, j);

    ASSERT_EQ(f(1), 16);
    ASSERT_EQ(f(-1), 14);
    ASSERT_EQ(f(5), 20);
}
