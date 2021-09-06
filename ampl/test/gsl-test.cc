/*
 Tests of the AMPL bindings for GNU Scientific Library.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */
#ifndef NOMINMAX
  #define NOMINMAX
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_version.h>

#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cstring>

#include "gtest/gtest.h"
// #define DEBUG_DIFFERENTIATOR
#include "function.h"
#include "asl.h"

using std::string;
using std::vector;

using fun::BitSet;
using fun::DERIVS;
using fun::DerivativeBinder;
using fun::Differentiator;
using fun::Function;
using fun::FunctionInfo;
using fun::HES;
using fun::MakeArgs;
using fun::Tuple;
using fun::Type;
using fun::Variant;

namespace fun {
template <>
gsl_sf_result *Convert<gsl_sf_result*>(const Variant &v) {
  return static_cast<gsl_sf_result*>(v.pointer());
}

template <>
const gsl_rng *Convert<const gsl_rng*>(const Variant &v) {
  return static_cast<const gsl_rng*>(v.pointer());
}
}

extern "C" void funcadd(AmplExports *ae);

namespace {

// Converts error estimate returned by Diff into an absolute tolerance to
// be used in EXPECT_NEAR.
double ConvertErrorToTolerance(double error) {
  return error != 0 ? error * 1000 : 1e-10;
}

class Error {
 private:
  string str_;

 public:
  explicit Error(const string &s) : str_(s) {}
  operator const char*() const { return str_.c_str(); }
  const char* c_str() const { return str_.c_str(); }
};

FunctionInfo::Result EvalError(const Function &f, const Tuple &args,
                               char prefix = 0, const char *suffix = "") {
  std::ostringstream os;
  if (prefix)
    os << prefix;
  os << "can't evaluate " << f.name() << suffix << args;
  return FunctionInfo::Result(os.str().c_str());
}

FunctionInfo::Result DerivError(const Function &f, const Tuple &args) {
  return EvalError(f, args, '\'', "'");
}

FunctionInfo::Result HesError(const Function &f, const Tuple &args) {
  return EvalError(f, args, '"', "''");
}

Error NotIntError(const string &arg_name, double value = 0.5) {
  std::ostringstream os;
  os << "argument '" << arg_name
      << "' can't be represented as int, " << arg_name << " = " << value;
  return Error(os.str());
}

Error NotUIntError(const string &arg_name, double value) {
  std::ostringstream os;
  os << "argument '" << arg_name
    << "' can't be represented as unsigned int, " << arg_name << " = " << value;
  return Error(os.str());
}

struct NoDeriv : FunctionInfo {
  explicit NoDeriv(const char *arg_names = "") : FunctionInfo(arg_names) {}

  Result GetDerivative(const Function &, unsigned, const Tuple &) const {
    return Result("'derivatives are not provided");
  }
};

#define EXPECT_ERROR(expected_message, result) \
  EXPECT_STREQ(expected_message, (result).error())

// Check if the value returned by af is correct.
void CheckFunction(double value, const Function &f, const Tuple &args) {
  std::ostringstream os;
  os << "Checking if " << f.name() << args << " = " << value;
  SCOPED_TRACE(os.str());
  if (gsl_isnan(value))
    EXPECT_ERROR(EvalError(f, args).error(), f(args));
  else
    EXPECT_EQ(value, f(args).value()) << f.name() << args;
}

typedef double (*FuncND)(int, double);

const double POINTS[] = {-5, -2, -1.23, -1, 0, 1, 1.23, 2, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

class GSLTest : public ::testing::Test {
 public:
  GSLTest() {}

 protected:
  static const FunctionInfo info;  // Default function info.

  Differentiator diff;

  // Differentiator statistics.
  struct Stats {
    unsigned num_calls_;
    unsigned num_nans_;
    double min_error_;
    double max_error_;

    Stats() :
      num_calls_(0), num_nans_(0),
      min_error_(std::numeric_limits<double>::max()),
      max_error_(std::numeric_limits<double>::min()) {}

    ~Stats() {
      std::cout << "Called numerical differentiation "
          << num_calls_ << " times with " << num_nans_ << " "
          << (num_nans_ == 1 ? "NaN" : "NaNs")
          << " detected" << std::endl;
      std::cout << "Error min=" << min_error_;
      std::cout << "  max=" << max_error_ << std::endl;
    }
  };
  static Stats stats_;

  template <typename F>
  double Diff(F f, double x, double *error = 0);

  static fun::Library lib_;
  static gsl_rng *rng_;

  static void SetUpTestCase() {
    lib_.Load();
    rng_ = gsl_rng_alloc(gsl_rng_default);

    // Check that the error handler is off.
    Function f(GetFunction("gsl_sf_airy_zero_Ai"));
    if (!f(-1).error())
      throw std::runtime_error("expected error");

    // Turn off the handler in the test in case it is not linked to the same
    // GSL library instance.
    gsl_set_error_handler_off();
  }

  static void TearDownTestCase() {
    gsl_rng_free(rng_);
  }

  // Returns an AMPL function by name.
  static Function GetFunction(const char *name, const FunctionInfo &info) {
    const func_info *fi = lib_.GetFunction(name);
    if (!fi)
      throw std::runtime_error(string("function not found: ") + name);
    return Function(&lib_, fi, &info);
  }

  static Function GetFunction(const char *name) {
    return GetFunction(name, info);
  }

  static bool CheckDerivative(
      double deriv, double numerical_deriv, double error);

  template <typename F>
  bool CheckDerivative(F f, const Function &af,
                       unsigned arg_index, const Tuple &args);

  static const unsigned NO_ARG = ~0u;

  template <typename F>
  void CheckSecondDerivatives(F f, const Function &af,
      const Tuple &args, unsigned skip_arg = NO_ARG);

  // Tests a function taking an integer and a double parameter.
  // test_x is a value of x where the function can be computed for very large
  // and very small n. If there is no such x or it is not known, then test_x
  // should be GSL_NAN.
  void TestFuncND(const Function &af, FuncND f, double test_x);

  template <typename F>
  void TestFunc(const Function &af, F f, Tuple &args, unsigned arg_index);

  template <typename F>
  void DoTestFunc(const Function &af, F f) {
    Tuple args(f.GetNumArgs());
    TestFunc(af, f, args, 0);
  }

  template <typename F>
  void TestFuncBindPrec(const Function &af, F f);

  template <typename F>
  void TestFuncBindPointers(const Function &af, F f);

  template <typename F>
  void TestFunc(const Function &af, F f) {
    TestFuncBindPointers(af, fun::FunctionPointer(f));
  }
};

GSLTest::Stats GSLTest::stats_;

const FunctionInfo GSLTest::info;
fun::Library GSLTest::lib_(AMPLGSL_DLL_NAME);
gsl_rng *GSLTest::rng_;

template <typename F>
double GSLTest::Diff(F f, double x, double *error) {
  ++stats_.num_calls_;
  double err = GSL_NAN;
  bool detected_nan = false;
  double deriv = diff(f, x, &err, &detected_nan);
  if (!gsl_isnan(deriv)) {
    stats_.min_error_ = std::min(stats_.min_error_, err);
    stats_.max_error_ = std::max(stats_.max_error_, err);
  }
  if (detected_nan)
    ++stats_.num_nans_;
  if (error)
    *error = err;
  return deriv;
}

bool GSLTest::CheckDerivative(
    double deriv, double numerical_deriv, double error) {
  if (deriv == numerical_deriv)
    return true;
  double abs_tolerance = ConvertErrorToTolerance(error);
  double diff = fabs(deriv - numerical_deriv);
  if (diff <= abs_tolerance)
    return true;
  if (fabs(numerical_deriv) > 0.001) {
    // If the derivative is not too close to zero, compare using a
    // relative tolerance.
    double relative_error = fabs(diff / numerical_deriv);
    if (relative_error < 1e-6)
      return false;
  }
  EXPECT_NEAR(numerical_deriv, deriv, abs_tolerance);
  return true;
}

// Checks if the value of the derivative returned by af agrees with the
// value returned by numerical differentiation of function f.
// arg_index: index of an argument with respect to which to differentiate
// args: point at which the derivative is computed
template <typename F>
bool GSLTest::CheckDerivative(F f, const Function &af,
                              unsigned arg_index, const Tuple &args) {
  std::ostringstream os;
  os << "Checking d/dx" << arg_index << " " << af.name() << " at " << args;
  SCOPED_TRACE(os.str());

  BitSet use_deriv(args.size(), false);
  use_deriv[arg_index] = true;

  double error = 0;
  double x = args[arg_index].number();
  FunctionInfo::Result deriv_result = af.GetDerivative(arg_index, args);
  double numerical_deriv = 0;
  if (!deriv_result.error()) {
    numerical_deriv = Diff(f, x, &error);
    double overridden_deriv = deriv_result.value();
    if (!gsl_isnan(overridden_deriv) && overridden_deriv != numerical_deriv) {
      std::cout << "Overriding d/dx" << arg_index << " " << af.name()
        << " at " << args << ", computed = " << numerical_deriv
        << ", overridden = " << overridden_deriv << std::endl;
      numerical_deriv = overridden_deriv;
    }
    if (!gsl_isnan(numerical_deriv)) {
      double deriv = af(args, DERIVS, use_deriv).deriv(arg_index);
      if (!CheckDerivative(deriv, numerical_deriv, error)) {
        std::cout << "Absolute tolerance of " << ConvertErrorToTolerance(error)
          << " not reached for d/dx" << arg_index << " " << af.name()
          << " at " << args << ", computed = " << numerical_deriv
          << ", actual = " << deriv << std::endl;
      }
      return true;
    }
  }

  // Call f and af equal number of times to keep RNGs in sync.
  double value = f(x);
  Function::Result r = af(args, DERIVS, use_deriv);
  if (!gsl_isnan(value)) {
    if (deriv_result.error())
      EXPECT_ERROR(deriv_result.error(), r);
    else
      EXPECT_ERROR(DerivError(af, args).error(), r);
  } else {
    EXPECT_TRUE(r.error() != 0);
  }
  return false;
}

// Checks if the values of the second partial derivatives returned by af
// agree with the values returned by numerical differentiation of the first
// partial derivatives.
// args: point at which the derivatives are computed
// skip_arg: index of an argument with respect to which not to differentiate
template <typename F>
void GSLTest::CheckSecondDerivatives(F f, const Function &af,
    const Tuple &args, unsigned skip_arg) {
  unsigned num_args = static_cast<unsigned>(args.size());
  if (skip_arg == NO_ARG) {
    for (unsigned i = 0; i < num_args; ++i) {
      if (af.GetDerivative(i, args).error()) {
        skip_arg = i;
        break;
      }
    }
  }
  for (unsigned i = 0; i < num_args; ++i) {
    if (i == skip_arg) continue;
    for (unsigned j = 0; j < num_args; ++j) {
      if (j == skip_arg) continue;
      BitSet use_deriv(num_args, false);
      use_deriv[i] = true;
      use_deriv[j] = true;
      double error = 0;
      FunctionInfo::Result deriv_result = af.GetSecondDerivative(i, j, args);
      double d = GSL_NAN;
      // Call f and af equal number of times to keep RNGs in sync.
      f(args);
      bool no_first_deriv = af(args, DERIVS, use_deriv).error() != 0;
      if (!deriv_result.error() && !no_first_deriv) {
        d = Diff(DerivativeBinder(af, j, i, args), args[i].number(), &error);
        double overridden_deriv = deriv_result.value();
        if (!gsl_isnan(overridden_deriv) && overridden_deriv != d) {
          std::cout << "Overriding d/dx" << i << " d/dx" << j << " "
            << af.name() << " at " << args << ", computed = " << d
            << ", overridden = " << overridden_deriv << std::endl;
          d = overridden_deriv;
        }
        if (!gsl_isnan(d)) {
          unsigned ii = i, jj = j;
          if (ii > jj) std::swap(ii, jj);
          unsigned hes_index = ii * (2 * num_args - ii - 1) / 2 + jj;
          double actual_deriv = af(args, HES, use_deriv).hes(hes_index);
          if (!CheckDerivative(actual_deriv, d, error)) {
            std::cout << "Absolute tolerance of "
              << ConvertErrorToTolerance(error)
              << " not reached for d/dx" << i << " d/dx" << j << " "
              << af.name() << " at " << args << ", computed = " << d
              << ", actual = " << actual_deriv << std::endl;
          }
          continue;
        }
      }
      // Call f and af equal number of times to keep RNGs in sync.
      f(args);
      Function::Result r = af(args, HES, use_deriv);
      if (no_first_deriv)
        EXPECT_TRUE(r.error() != 0);
      else if (deriv_result.error())
        EXPECT_ERROR(deriv_result.error(), r);
      else
        EXPECT_ERROR(HesError(af, args).error(), r);
    }
  }
}

void GSLTest::TestFuncND(const Function &af, FuncND f, double test_x) {
  TestFunc(af, f);

  // These points are tested separately because of various problems, e.g.
  // gsl_sf_bessel_Jn(n, x) and gsl_sf_bessel_In(n, x) take too much time
  // (hang?) for n = INT_MIN and gsl_sf_bessel_Yn(n, x) returns different
  // values close to x = 0 when called different times for n = INT_MIN.

  BitSet use_deriv("01");
  string arg_name(af.GetArgName(0));
  EXPECT_ERROR(
      ("'can't compute derivative: argument '" + arg_name + "' too small, " +
      arg_name + " = -2147483648").c_str(),
      af(MakeArgs(INT_MIN, test_x), DERIVS, use_deriv));
  EXPECT_ERROR(
      ("'can't compute derivative: argument '" + arg_name + "' too large, " +
      arg_name + " = 2147483647").c_str(),
      af(MakeArgs(INT_MAX, test_x), DERIVS, use_deriv));
  if (!af(MakeArgs(INT_MIN + 1, test_x)).error()) {
    EXPECT_TRUE(!gsl_isnan(af(MakeArgs(INT_MIN + 1, test_x), DERIVS,
        use_deriv).deriv(1)));
  }
  if (!af(MakeArgs(INT_MAX - 1, test_x)).error()) {
    EXPECT_TRUE(!gsl_isnan(af(MakeArgs(INT_MAX - 1, test_x), DERIVS,
        use_deriv).deriv(1)));
  }

  EXPECT_ERROR(
      ("'can't compute derivative: argument '" + arg_name + "' too small, " +
      arg_name + " = -2147483647").c_str(),
      af(MakeArgs(INT_MIN + 1, test_x), HES, use_deriv));
  EXPECT_ERROR(
      ("'can't compute derivative: argument '" + arg_name + "' too large, " +
      arg_name + " = 2147483646").c_str(),
      af(MakeArgs(INT_MAX - 1, test_x), HES, use_deriv));
  if (!af(MakeArgs(INT_MIN + 2, test_x)).error()) {
    EXPECT_TRUE(!gsl_isnan(af(MakeArgs(INT_MIN + 2, test_x), HES,
        use_deriv).hes(2)));
  }
  if (!af(MakeArgs(INT_MAX - 2, test_x)).error()) {
    EXPECT_TRUE(!gsl_isnan(af(MakeArgs(INT_MAX - 2, test_x), HES,
        use_deriv).hes(2)));
  }
}

// Checks that the error is returned when trying to get a derivative
// with respect to an integer argument.
// Returns true if the argument is integer, false otherwise.
template <typename F>
bool CheckArg(F f, const Function &af, const Tuple &args, unsigned arg_index) {
  if (f.GetArgType(arg_index) == fun::DOUBLE)
    return false;
  std::ostringstream os;
  BitSet use_deriv(args.size(), false);
  use_deriv[arg_index] = true;
  os << "'argument '" << af.GetArgName(arg_index) << "' is not constant";
  // Call f and af equal number of times to keep RNGs in sync.
  f(args);
  EXPECT_STREQ(os.str().c_str(), af(args, DERIVS, use_deriv).error());
  return true;
}

template <typename F>
void GSLTest::TestFunc(
    const Function &af, F f, Tuple &args, unsigned arg_index) {
  unsigned num_args = static_cast<unsigned>(args.size());
  if (arg_index == 0) {
    bool has_double_arg = false;
    for (unsigned i = 0; i < num_args; ++i) {
      if (f.GetArgType(i) == fun::DOUBLE) {
        has_double_arg = true;
        args[i] = Variant::FromDouble(GSL_NAN);
      } else {
        args[i] = Variant::FromDouble(0);
      }
    }
    if (has_double_arg)
      EXPECT_ERROR(EvalError(af, args).error(), af(args));
  }
  if (arg_index < num_args) {
    for (size_t i = 0; i != NUM_POINTS; ++i) {
      args[arg_index] = Variant::FromDouble(POINTS[i]);
      TestFunc(af, f, args, arg_index + 1);
    }
    return;
  }
  for (unsigned i = 0; i < num_args; ++i) {
    double arg = args[i].number();
    if (f.GetArgType(i) == fun::INT && static_cast<int>(arg) != arg) {
      EXPECT_STREQ(NotIntError(af.GetArgName(i), arg), af(args).error());
      return;
    }
    if (f.GetArgType(i) == fun::UINT &&
        static_cast<unsigned>(arg) != arg) {
      EXPECT_STREQ(NotUIntError(af.GetArgName(i), arg), af(args).error());
      return;
    }
  }
  if (af.SkipPoint(args))
    return;
  CheckFunction(f(args), af, args);
  for (unsigned i = 0; i < num_args; ++i) {
    // Call f and af equal number of times to keep RNGs in sync.
    BitSet use_deriv(num_args, false);
    Function::Result result = af(args, DERIVS, use_deriv);
    use_deriv[i] = true;
    if (gsl_isnan(f(args)))
      EXPECT_ERROR(EvalError(af, args).error(), result);
    else if (!CheckArg(f, af, args, i))
      CheckDerivative(BindAllButOne(f, args, i), af, i, args);
  }
  CheckSecondDerivatives(f, af, args);
}

template <typename F>
class ResultBinder {
 private:
  F f_;

 public:
  explicit ResultBinder(F f) : f_(f) {}

  unsigned GetNumArgs() const { return f_.GetNumArgs() - 1; }

  Type GetArgType(unsigned arg_index) const {
    if (arg_index >= GetNumArgs())
      throw std::out_of_range("argument index is out of range");
    return f_.GetArgType(arg_index);
  }

  double operator()(const Tuple &args) const {
    gsl_sf_result result = gsl_sf_result();
    return fun::OneBinder<F, gsl_sf_result*, int>(
        f_, &result, GetNumArgs())(args) ? GSL_NAN : result.val;
  }
};

template <typename F>
ResultBinder<F> BindResult(F f) {
  return ResultBinder<F>(f);
}

template <typename F>
void GSLTest::TestFuncBindPrec(const Function &af, F f) {
  unsigned num_args = f.GetNumArgs();
  if (af.nargs() == static_cast<int>(num_args) - 1 &&
      f.GetArgType(num_args - 1) == fun::UINT) {
    // If the last argument is a mode, bind it to GSL_PREC_DOUBLE.
    DoTestFunc(af, BindOne(f, GSL_PREC_DOUBLE, num_args - 1));
    return;
  }
  DoTestFunc(af, f);
}

template <typename F>
void GSLTest::TestFuncBindPointers(const Function &af, F f) {
  unsigned num_args = f.GetNumArgs();
  if (af.nargs() < static_cast<int>(num_args)) {
    fun::Type first_type = f.GetArgType(0);
    if (af.ftype() == FUNCADD_RANDOM_VALUED && first_type == fun::POINTER) {
      // If the first argument is a random number generator, bind it.
      TestFuncBindPrec(af, BindOne(f, rng_, 0));
      return;
    }
    if (f.GetArgType(num_args - 1) == fun::POINTER) {
      // If the last argument is a result, bind it.
      TestFuncBindPrec(af, BindResult(f));
      return;
    }
  }
  TestFuncBindPrec(af, f);
}

#define TEST_FUNC(name) TestFunc(GetFunction(#name, this->info), name)
#define TEST_FUNC2(name, info) TestFunc(GetFunction(#name, info), name)

#define TEST_EFUNC(name) TestFunc(GetFunction(#name, this->info), name##_e)
#define TEST_EFUNC2(name, info) TestFunc(GetFunction(#name, info), name##_e)

#define TEST_FUNC_ND(name, test_x) \
  TestFuncND(GetFunction(#name, FunctionInfo("n")), name, test_x)

TEST_F(GSLTest, NoLibErrors) {
  EXPECT_EQ("", lib_.error());
}

TEST_F(GSLTest, Elementary) {
  TEST_FUNC(gsl_log1p);
  TEST_FUNC(gsl_expm1);
  TEST_FUNC(gsl_hypot);
  TEST_FUNC(gsl_hypot3);
}

TEST_F(GSLTest, AiryA) {
  TEST_EFUNC(gsl_sf_airy_Ai);
  TEST_EFUNC(gsl_sf_airy_Ai_scaled);
}

TEST_F(GSLTest, AiryB) {
  TEST_EFUNC(gsl_sf_airy_Bi);
  TEST_EFUNC(gsl_sf_airy_Bi_scaled);
}

TEST_F(GSLTest, AiryZero) {
  FunctionInfo info("s");
  TEST_EFUNC2(gsl_sf_airy_zero_Ai, info);
  TEST_EFUNC2(gsl_sf_airy_zero_Bi, info);
  TEST_EFUNC2(gsl_sf_airy_zero_Ai_deriv, info);
  TEST_EFUNC2(gsl_sf_airy_zero_Bi_deriv, info);
}

TEST_F(GSLTest, BesselJ) {
  TEST_EFUNC(gsl_sf_bessel_J0);
  TEST_EFUNC(gsl_sf_bessel_J1);
  TEST_FUNC_ND(gsl_sf_bessel_Jn, 0);
}

TEST_F(GSLTest, BesselY) {
  TEST_EFUNC(gsl_sf_bessel_Y0);
  TEST_EFUNC(gsl_sf_bessel_Y1);
  TEST_FUNC_ND(gsl_sf_bessel_Yn, 1);
}

TEST_F(GSLTest, BesselI) {
  TEST_EFUNC(gsl_sf_bessel_I0);
  TEST_EFUNC(gsl_sf_bessel_I1);
  TEST_FUNC_ND(gsl_sf_bessel_In, 0);
  TEST_EFUNC(gsl_sf_bessel_I0_scaled);
  TEST_EFUNC(gsl_sf_bessel_I1_scaled);
  TEST_FUNC_ND(gsl_sf_bessel_In_scaled, 0);
}

TEST_F(GSLTest, BesselK) {
  TEST_EFUNC(gsl_sf_bessel_K0);
  TEST_EFUNC(gsl_sf_bessel_K1);
  TEST_FUNC_ND(gsl_sf_bessel_Kn, 1);
  TEST_EFUNC(gsl_sf_bessel_K0_scaled);
  TEST_EFUNC(gsl_sf_bessel_K1_scaled);
  TEST_FUNC_ND(gsl_sf_bessel_Kn_scaled, 1);
}

TEST_F(GSLTest, Besselj) {
  TEST_EFUNC(gsl_sf_bessel_j0);
  TEST_EFUNC(gsl_sf_bessel_j1);
  TEST_EFUNC(gsl_sf_bessel_j2);
  TEST_EFUNC2(gsl_sf_bessel_jl, FunctionInfo("l x"));
}

TEST_F(GSLTest, Bessely) {
  TEST_EFUNC(gsl_sf_bessel_y0);
  TEST_EFUNC(gsl_sf_bessel_y1);
  TEST_EFUNC(gsl_sf_bessel_y2);
  TEST_EFUNC2(gsl_sf_bessel_yl, FunctionInfo("l x"));
}

TEST_F(GSLTest, Besseli) {
  TEST_EFUNC(gsl_sf_bessel_i0_scaled);
  TEST_EFUNC(gsl_sf_bessel_i1_scaled);
  TEST_EFUNC(gsl_sf_bessel_i2_scaled);
  TEST_EFUNC2(gsl_sf_bessel_il_scaled, FunctionInfo("l x"));
}

TEST_F(GSLTest, Besselk) {
  TEST_EFUNC(gsl_sf_bessel_k0_scaled);
  TEST_EFUNC(gsl_sf_bessel_k1_scaled);
  TEST_EFUNC(gsl_sf_bessel_k2_scaled);
  TEST_EFUNC2(gsl_sf_bessel_kl_scaled, FunctionInfo("l x"));
}

struct BesselFractionalOrderInfo : FunctionInfo {
  Result GetDerivative(const Function &f,
      unsigned arg_index, const Tuple &args) const {
    // Partial derivative with respect to nu is not provided.
    if (arg_index == 0)
      return Result("'argument 'nu' is not constant");
    // Computing gsl_sf_bessel_*nu'(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 1, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 1.
    if (args[0].number() < 1)
      return DerivError(f, args);
    return Result();
  }

  Result GetSecondDerivative(
      const Function &f, unsigned, unsigned, const Tuple &args) const {
    // Computing gsl_sf_bessel_*nu''(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 2, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 2.
    return args[0].number() < 2 ? HesError(f, args) : Result();
  }
};

TEST_F(GSLTest, BesselFractionalOrder) {
  BesselFractionalOrderInfo info;
  // TEST_EFUNC2(gsl_sf_bessel_Jnu, info); // FIXME: (GSL v2.5.0) fails with arguments from POINTS
  // TEST_EFUNC2(gsl_sf_bessel_Ynu, info); // FIXME: (GSL v2.5.0) fails with arguments from POINTS
  TEST_EFUNC2(gsl_sf_bessel_Inu, info);
  TEST_EFUNC2(gsl_sf_bessel_Inu_scaled, info);
  TEST_EFUNC2(gsl_sf_bessel_Knu, info);
  TEST_EFUNC2(gsl_sf_bessel_lnKnu, info);
  TEST_EFUNC2(gsl_sf_bessel_Knu_scaled, info);
}

TEST_F(GSLTest, BesselZero) {
  TEST_EFUNC2(gsl_sf_bessel_zero_J0, NoDeriv("s"));
  TEST_EFUNC2(gsl_sf_bessel_zero_J1, NoDeriv("s"));
  TEST_EFUNC2(gsl_sf_bessel_zero_Jnu, NoDeriv("nu s"));
}

struct ClausenFunctionInfo : FunctionInfo {
  Result GetDerivative(const Function &, unsigned, const Tuple &args) const {
    return Result(args[0].number() == 0 ? GSL_POSINF : GSL_NAN);
  }
};

TEST_F(GSLTest, Clausen) {
  TEST_EFUNC2(gsl_sf_clausen, ClausenFunctionInfo());
}

TEST_F(GSLTest, Hydrogenic) {
  TEST_EFUNC(gsl_sf_hydrogenicR_1);
  TEST_EFUNC2(gsl_sf_hydrogenicR, NoDeriv("n l z r"));
}

TEST_F(GSLTest, Coulomb) {
  TEST_EFUNC2(gsl_sf_coulomb_CL, NoDeriv());
}

TEST_F(GSLTest, Coupling3j) {
  double value = gsl_sf_coupling_3j(8, 20, 12, -2, 12, -10);
  EXPECT_NEAR(0.0812695955, value, 1e-5);
  Function f = GetFunction("gsl_sf_coupling_3j");
  Tuple args(MakeArgs(8, 20, 12, -2, 12, -10));
  EXPECT_EQ(value, f(args).value());
  f(MakeArgs(0, 0, 0, 0, 0, 0));
  EXPECT_ERROR(NotIntError("two_ja"), f(MakeArgs(0.5, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jb"), f(MakeArgs(0, 0.5, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jc"), f(MakeArgs(0, 0, 0.5, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_ma"), f(MakeArgs(0, 0, 0, 0.5, 0, 0)));
  EXPECT_ERROR(NotIntError("two_mb"), f(MakeArgs(0, 0, 0, 0, 0.5, 0)));
  EXPECT_ERROR(NotIntError("two_mc"), f(MakeArgs(0, 0, 0, 0, 0, 0.5)));
  const char *error = "'argument 'two_ja' is not constant";
  EXPECT_ERROR(error, f(args, DERIVS));
  EXPECT_ERROR(error, f(args, HES));
}

TEST_F(GSLTest, Coupling6j) {
  double value = gsl_sf_coupling_6j(2, 4, 6, 8, 10, 12);
  EXPECT_NEAR(0.0176295295, value, 1e-7);
  Function f = GetFunction("gsl_sf_coupling_6j");
  Tuple args(MakeArgs(2, 4, 6, 8, 10, 12));
  EXPECT_EQ(value, f(args).value());
  EXPECT_TRUE(f(MakeArgs(0, 0, 0, 0, 0, 0)).error() == 0);
  EXPECT_ERROR(NotIntError("two_ja"), f(MakeArgs(0.5, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jb"), f(MakeArgs(0, 0.5, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jc"), f(MakeArgs(0, 0, 0.5, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jd"), f(MakeArgs(0, 0, 0, 0.5, 0, 0)));
  EXPECT_ERROR(NotIntError("two_je"), f(MakeArgs(0, 0, 0, 0, 0.5, 0)));
  EXPECT_ERROR(NotIntError("two_jf"), f(MakeArgs(0, 0, 0, 0, 0, 0.5)));
  const char *error = "'argument 'two_ja' is not constant";
  EXPECT_ERROR(error, f(args, DERIVS));
  EXPECT_ERROR(error, f(args, HES));
}

TEST_F(GSLTest, Coupling9j) {
  double value = gsl_sf_coupling_9j(6, 16, 18, 8, 20, 14, 12, 10, 4);
  EXPECT_NEAR(-0.000775648399, value, 1e-9);
  Function f = GetFunction("gsl_sf_coupling_9j");
  Tuple args(MakeArgs(6, 16, 18, 8, 20, 14, 12, 10, 4));
  EXPECT_EQ(value, f(args).value());
  EXPECT_TRUE(f(MakeArgs(0, 0, 0, 0, 0, 0, 0, 0, 0)).error() == 0);
  EXPECT_ERROR(NotIntError("two_ja"), f(MakeArgs(0.5, 0, 0, 0, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jb"), f(MakeArgs(0, 0.5, 0, 0, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jc"), f(MakeArgs(0, 0, 0.5, 0, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jd"), f(MakeArgs(0, 0, 0, 0.5, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_je"), f(MakeArgs(0, 0, 0, 0, 0.5, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jf"), f(MakeArgs(0, 0, 0, 0, 0, 0.5, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jg"), f(MakeArgs(0, 0, 0, 0, 0, 0, 0.5, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jh"), f(MakeArgs(0, 0, 0, 0, 0, 0, 0, 0.5, 0)));
  EXPECT_ERROR(NotIntError("two_ji"), f(MakeArgs(0, 0, 0, 0, 0, 0, 0, 0, 0.5)));
  const char *error = "'argument 'two_ja' is not constant";
  EXPECT_ERROR(error, f(args, DERIVS));
  EXPECT_ERROR(error, f(args, HES));
}

TEST_F(GSLTest, Dawson) {
  TEST_EFUNC(gsl_sf_dawson);
}

TEST_F(GSLTest, Debye) {
  TEST_EFUNC(gsl_sf_debye_1);
  TEST_EFUNC(gsl_sf_debye_2);
  TEST_EFUNC(gsl_sf_debye_3);
  TEST_EFUNC(gsl_sf_debye_4);
  TEST_EFUNC(gsl_sf_debye_5);
  TEST_EFUNC(gsl_sf_debye_6);
}

struct DilogFunctionInfo : FunctionInfo {
  Result GetDerivative(const Function &, unsigned, const Tuple &args) const {
    return Result(args[0].number() == 1 ? GSL_POSINF : GSL_NAN);
  }
};

TEST_F(GSLTest, Dilog) {
  TEST_EFUNC2(gsl_sf_dilog, DilogFunctionInfo());
}

int ellint_D(double phi, double k, gsl_mode_t mode, gsl_sf_result *result) {
#if GSL_MAJOR_VERSION >= 2
  return gsl_sf_ellint_D_e(phi, k, mode, result);
#else
  return gsl_sf_ellint_D_e(phi, k, 0, mode, result);
#endif
}

TEST_F(GSLTest, EllInt) {
  TEST_EFUNC(gsl_sf_ellint_Kcomp);
  TEST_EFUNC(gsl_sf_ellint_Ecomp);
  TEST_EFUNC(gsl_sf_ellint_Pcomp);
  TEST_EFUNC(gsl_sf_ellint_F);
  TEST_EFUNC(gsl_sf_ellint_E);
  TEST_EFUNC2(gsl_sf_ellint_P, NoDeriv());
  TestFunc(GetFunction("gsl_sf_ellint_D", NoDeriv()), ellint_D);
  TEST_EFUNC2(gsl_sf_ellint_RC, NoDeriv());
  TEST_EFUNC2(gsl_sf_ellint_RD, NoDeriv());
  TEST_EFUNC2(gsl_sf_ellint_RF, NoDeriv());
  TEST_EFUNC2(gsl_sf_ellint_RJ, NoDeriv());
}

TEST_F(GSLTest, Erf) {
  TEST_EFUNC(gsl_sf_erf);
  TEST_EFUNC(gsl_sf_erfc);
  TEST_EFUNC(gsl_sf_log_erfc);
  TEST_EFUNC(gsl_sf_erf_Z);
  TEST_EFUNC(gsl_sf_erf_Q);
  TEST_EFUNC(gsl_sf_hazard);
}

TEST_F(GSLTest, ExpInt) {
  TEST_EFUNC(gsl_sf_expint_E1);
  TEST_EFUNC(gsl_sf_expint_E2);
  TEST_EFUNC2(gsl_sf_expint_En, FunctionInfo("n x"));
  TEST_EFUNC(gsl_sf_expint_Ei);
  TEST_EFUNC(gsl_sf_Shi);
  TEST_EFUNC(gsl_sf_Chi);
  TEST_EFUNC(gsl_sf_expint_3);
  TEST_EFUNC(gsl_sf_Si);
  TEST_EFUNC(gsl_sf_Ci);
  TEST_EFUNC(gsl_sf_atanint);
}

TEST_F(GSLTest, FermiDirac) {
  TEST_EFUNC(gsl_sf_fermi_dirac_m1);
  TEST_EFUNC(gsl_sf_fermi_dirac_0);
  TEST_EFUNC(gsl_sf_fermi_dirac_1);
  TEST_EFUNC(gsl_sf_fermi_dirac_2);
  TEST_EFUNC2(gsl_sf_fermi_dirac_int, FunctionInfo("j x"));
  TEST_EFUNC2(gsl_sf_fermi_dirac_mhalf, NoDeriv());
  TEST_EFUNC2(gsl_sf_fermi_dirac_half, NoDeriv());
  TEST_EFUNC(gsl_sf_fermi_dirac_3half);
  TEST_EFUNC(gsl_sf_fermi_dirac_inc_0);
}

struct LnGammaInfo : FunctionInfo {
  Result GetDerivative(const Function &af, unsigned , const Tuple &args) const {
    double x = args[0].number();
    return x == -1 || x == -2 ? DerivError(af, args) : Result();
  }
};

TEST_F(GSLTest, Gamma) {
  Function gamma = GetFunction("gsl_sf_gamma");
  EXPECT_NEAR(-0.129354, gamma(-0.5, DERIVS).deriv(), 1e-6);
  EXPECT_NEAR(-31.6778, gamma(-0.5, HES).hes(), 1e-4);
  EXPECT_NEAR(1.19786e100, gamma(71).value(), 1e95);
  EXPECT_NE(0, gsl_isinf(gamma(1000).value()));
  TEST_EFUNC(gsl_sf_gamma);
  TEST_EFUNC2(gsl_sf_lngamma, LnGammaInfo());
  TEST_EFUNC(gsl_sf_gammastar);
  TEST_EFUNC(gsl_sf_gammainv);
}

TEST_F(GSLTest, Poch) {
  TEST_EFUNC2(gsl_sf_poch, NoDeriv());
  TEST_EFUNC2(gsl_sf_lnpoch, NoDeriv());
  TEST_EFUNC2(gsl_sf_pochrel, NoDeriv());
}

struct GammaIncInfo : FunctionInfo {
  Result GetDerivative(
      const Function &af, unsigned arg_index, const Tuple &args) const {
    // Partial derivative with respect to a is not provided.
    if (arg_index == 0)
      return Result("'argument 'a' is not constant");
    if (args[1].number() == 0)
      return DerivError(af, args);
    return Result();
  }
};

TEST_F(GSLTest, GammaInc) {
  TEST_EFUNC2(gsl_sf_gamma_inc, GammaIncInfo());
  TEST_EFUNC2(gsl_sf_gamma_inc_Q, NoDeriv());
  TEST_EFUNC2(gsl_sf_gamma_inc_P, NoDeriv());
}

struct BetaInfo : FunctionInfo {
  Result GetDerivative(
      const Function &f, unsigned arg_index, const Tuple &args) const {
    if (gsl_isnan(gsl_sf_psi(args[0].number() + args[1].number())) ||
        args[arg_index].number() == 0) {
      return DerivError(f, args);
    }
    return Result();
  }
};

TEST_F(GSLTest, Beta) {
  TEST_EFUNC2(gsl_sf_beta, BetaInfo());
  TEST_EFUNC(gsl_sf_lnbeta);
  TEST_EFUNC2(gsl_sf_beta_inc, NoDeriv());
}

TEST_F(GSLTest, GegenPoly) {
  TEST_EFUNC(gsl_sf_gegenpoly_1);
  TEST_EFUNC(gsl_sf_gegenpoly_2);
  TEST_EFUNC(gsl_sf_gegenpoly_3);
  TEST_EFUNC2(gsl_sf_gegenpoly_n, NoDeriv("n"));
}

class NoDerivArg : public FunctionInfo {
 private:
  unsigned arg_index_;

 public:
  NoDerivArg(const char *arg_names, unsigned arg_index = 0)
  : FunctionInfo(arg_names), arg_index_(arg_index) {}

  Result GetDerivative(
      const Function &f, unsigned arg_index, const Tuple &) const {
    // Partial derivative with respect to the first argument is not provided.
    string error(string("'argument '") + f.GetArgName(arg_index_) +
        "' is not constant");
    return Result(arg_index == arg_index_ ? error.c_str() : "");
  }
};

struct Hyperg1F1Info : FunctionInfo {
  Result GetDerivative(const Function &f, unsigned, const Tuple &args) const {
    return args[1].number() <= 0 ? DerivError(f, args) : Result();
  }
};

TEST_F(GSLTest, Hyperg) {
  TEST_EFUNC2(gsl_sf_hyperg_0F1, NoDerivArg("c"));
  TEST_EFUNC2(gsl_sf_hyperg_1F1_int, Hyperg1F1Info().SetArgNames("m n x"));
  TEST_EFUNC2(gsl_sf_hyperg_1F1, NoDeriv());
  TEST_EFUNC2(gsl_sf_hyperg_U_int, NoDeriv("m n x"));
  TEST_EFUNC2(gsl_sf_hyperg_U, NoDeriv());
  TEST_EFUNC2(gsl_sf_hyperg_2F1, NoDeriv());
  TEST_EFUNC2(gsl_sf_hyperg_2F1_conj, NoDeriv());
  TEST_EFUNC2(gsl_sf_hyperg_2F1_renorm, NoDeriv());
  TEST_EFUNC2(gsl_sf_hyperg_2F1_conj_renorm, NoDeriv());
  TEST_EFUNC2(gsl_sf_hyperg_2F0, NoDeriv());
}

TEST_F(GSLTest, Laguerre) {
  TEST_EFUNC(gsl_sf_laguerre_1);
  TEST_EFUNC(gsl_sf_laguerre_2);
  TEST_EFUNC(gsl_sf_laguerre_3);
  TEST_EFUNC2(gsl_sf_laguerre_n, NoDeriv("n"));
}

struct LambertW0Info : FunctionInfo {
  Result GetDerivative(const Function &f, unsigned, const Tuple &args) const {
    return args[0].number() < -1 / M_E ? DerivError(f, args) : Result();
  }
};

struct LambertWm1Info : FunctionInfo {
  Result GetDerivative(const Function &f, unsigned, const Tuple &args) const {
    double x = args[0].number();
    return x < -1 / M_E || x == 0 ? DerivError(f, args) : Result();
  }
};

TEST_F(GSLTest, Lambert) {
  TEST_EFUNC2(gsl_sf_lambert_W0, LambertW0Info());
  TEST_EFUNC2(gsl_sf_lambert_Wm1, LambertWm1Info());
  EXPECT_NEAR(-13.8803,
      GetFunction("gsl_sf_lambert_Wm1")(-0.1, DERIVS).deriv(0), 1e-4);
}

TEST_F(GSLTest, Legendre) {
  TEST_EFUNC(gsl_sf_legendre_P1);
  TEST_EFUNC(gsl_sf_legendre_P2);
  TEST_EFUNC(gsl_sf_legendre_P3);
  TEST_EFUNC2(gsl_sf_legendre_Pl, FunctionInfo("l"));
  TEST_EFUNC(gsl_sf_legendre_Q0);
  TEST_EFUNC(gsl_sf_legendre_Q1);
  TEST_EFUNC2(gsl_sf_legendre_Ql, FunctionInfo("l"));
  TEST_EFUNC2(gsl_sf_legendre_Plm, NoDeriv("l m"));
  TEST_EFUNC2(gsl_sf_legendre_sphPlm, NoDeriv("l m"));
}

TEST_F(GSLTest, Conical) {
  TEST_EFUNC2(gsl_sf_conicalP_half, NoDeriv());
  TEST_EFUNC2(gsl_sf_conicalP_mhalf, NoDeriv());
  TEST_EFUNC2(gsl_sf_conicalP_0, NoDeriv());
  TEST_EFUNC2(gsl_sf_conicalP_1, NoDeriv());
  TEST_EFUNC2(gsl_sf_conicalP_sph_reg, NoDeriv("m"));
  TEST_EFUNC2(gsl_sf_conicalP_cyl_reg, NoDeriv("m"));
}

TEST_F(GSLTest, RadialLegendre) {
  TEST_EFUNC2(gsl_sf_legendre_H3d_0, NoDeriv());
  TEST_EFUNC2(gsl_sf_legendre_H3d_1, NoDeriv());
  TEST_EFUNC2(gsl_sf_legendre_H3d, NoDeriv("l"));
}

TEST_F(GSLTest, Log) {
  TEST_EFUNC(gsl_sf_log);
  TEST_EFUNC(gsl_sf_log_abs);
  TEST_EFUNC(gsl_sf_log_1plusx);
  TEST_EFUNC(gsl_sf_log_1plusx_mx);
}

TEST_F(GSLTest, Mathieu) {
#if GSL_MAJOR_VERSION >= 2
# define TEST_MATHIEU TEST_EFUNC2
#else
# define TEST_MATHIEU TEST_FUNC2
#endif
  TEST_MATHIEU(gsl_sf_mathieu_a, NoDeriv("n"));
  TEST_MATHIEU(gsl_sf_mathieu_b, NoDeriv("n"));
  TEST_MATHIEU(gsl_sf_mathieu_ce, NoDeriv("n"));
  TEST_MATHIEU(gsl_sf_mathieu_se, NoDeriv("n"));
  TEST_MATHIEU(gsl_sf_mathieu_Mc, NoDeriv("j n"));
  TEST_MATHIEU(gsl_sf_mathieu_Ms, NoDeriv("j n"));
}

TEST_F(GSLTest, Power) {
  TEST_EFUNC2(gsl_sf_pow_int, FunctionInfo("x n"));
}

struct DigammaInfo : FunctionInfo {
  Result GetDerivative(const Function &f, unsigned , const Tuple &args) const {
    double x = args[0].number();
    return x < 0 && ceil(x) == x ? DerivError(f, args) : Result();
  }
  Result GetSecondDerivative(
      const Function &f, unsigned , unsigned , const Tuple &args) const {
    // gsl_sf_psi_n doesn't support negative x.
    return args[0].number() < 0 ? HesError(f, args) : Result();
  }
};

struct TrigammaInfo : FunctionInfo {
  Result GetDerivative(const Function &f, unsigned , const Tuple &args) const {
    // gsl_sf_psi_n doesn't support negative x.
    return args[0].number() < 0 ? DerivError(f, args) : Result();
  }
};

struct PolygammaInfo : FunctionInfo {
  Result GetDerivative(
      const Function &f, unsigned , const Tuple &args) const {
    double x = args[1].number();
    if (args[0].number() == 0)
      return x < 0 && ceil(x) == x ? DerivError(f, args) : Result();
    return x < 0 ? DerivError(f, args) : Result();
  }
  Result GetSecondDerivative(const Function &f,
      unsigned , unsigned , const Tuple &args) const {
    // gsl_sf_psi_n doesn't support negative x.
    return args[1].number() < 0 ? HesError(f, args) : Result();
  }
};

TEST_F(GSLTest, Digamma) {
  TEST_EFUNC2(gsl_sf_psi_int, FunctionInfo("n"));
  TEST_EFUNC2(gsl_sf_psi, DigammaInfo());
  TEST_EFUNC2(gsl_sf_psi_1piy, NoDeriv());
  TEST_EFUNC2(gsl_sf_psi_1_int, FunctionInfo("n"));
  TEST_EFUNC2(gsl_sf_psi_1, TrigammaInfo());
  TEST_EFUNC2(gsl_sf_psi_n, PolygammaInfo().SetArgNames("n"));
}

TEST_F(GSLTest, Synchrotron) {
  TEST_EFUNC2(gsl_sf_synchrotron_1, NoDeriv());
  TEST_EFUNC2(gsl_sf_synchrotron_2, NoDeriv());
}

TEST_F(GSLTest, Transport) {
  TEST_EFUNC(gsl_sf_transport_2);
  TEST_EFUNC(gsl_sf_transport_3);
  TEST_EFUNC(gsl_sf_transport_4);
  TEST_EFUNC(gsl_sf_transport_5);
}

TEST_F(GSLTest, Zeta) {
  TEST_EFUNC2(gsl_sf_zeta_int, FunctionInfo("n"));
  TEST_EFUNC2(gsl_sf_zeta, NoDeriv());
  TEST_EFUNC2(gsl_sf_zetam1_int, FunctionInfo("n"));
  TEST_EFUNC2(gsl_sf_zetam1, NoDeriv());
  TEST_EFUNC2(gsl_sf_hzeta, NoDeriv());
  TEST_EFUNC2(gsl_sf_eta_int, FunctionInfo("n"));
  TEST_EFUNC2(gsl_sf_eta, NoDeriv());
}

TEST_F(GSLTest, Gaussian) {
  TEST_FUNC2(gsl_ran_gaussian, NoDeriv());
  TEST_FUNC(gsl_ran_gaussian_pdf);
  TEST_FUNC2(gsl_ran_gaussian_ziggurat, NoDeriv());
  TEST_FUNC2(gsl_ran_gaussian_ratio_method, NoDeriv());
  TEST_FUNC2(gsl_ran_ugaussian, NoDeriv());
  TEST_FUNC(gsl_ran_ugaussian_pdf);
  TEST_FUNC2(gsl_ran_ugaussian_ratio_method, NoDeriv());

  struct GaussianPInfo : FunctionInfo {
    Result GetDerivative(
        const Function &f,unsigned , const Tuple &args) const {
      return args[1].number() == 0 ? DerivError(f, args) : Result();
    }
  };
  TEST_FUNC2(gsl_cdf_gaussian_P, GaussianPInfo());

  TEST_FUNC2(gsl_cdf_gaussian_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_gaussian_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_gaussian_Qinv, NoDeriv());
  TEST_FUNC(gsl_cdf_ugaussian_P);
  TEST_FUNC2(gsl_cdf_ugaussian_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_ugaussian_Qinv, NoDeriv());
}

TEST_F(GSLTest, UGaussianPInv) {
  struct UGaussianPInvInfo : FunctionInfo {
    Result GetDerivative(const Function &, unsigned , const Tuple &args) const {
      double x = args[0].number();
      return x == 0 || x == 1 ? Result(GSL_POSINF) : Result();
    }
    Result GetSecondDerivative(
        const Function &, unsigned, unsigned, const Tuple &args) const {
      double x = args[0].number();
      if (x == 0) return Result(GSL_NEGINF);
      return x == 1 ? Result(GSL_POSINF) : Result();
    }
  };
  TEST_FUNC2(gsl_cdf_ugaussian_Pinv, UGaussianPInvInfo());
  Function pinv = GetFunction("gsl_cdf_ugaussian_Pinv");
  EXPECT_NEAR(2.876103, pinv(0.3, DERIVS).deriv(), 1e-5);
  EXPECT_NEAR(-4.33782, pinv(0.3, HES).hes(), 1e-5);
}

TEST_F(GSLTest, GaussianTail) {
  TEST_FUNC2(gsl_ran_gaussian_tail, NoDeriv());
  TEST_FUNC2(gsl_ran_gaussian_tail_pdf, NoDeriv());
  TEST_FUNC2(gsl_ran_ugaussian_tail, NoDeriv());
  TEST_FUNC2(gsl_ran_ugaussian_tail_pdf, NoDeriv());
}

TEST_F(GSLTest, Exponential) {
  TEST_FUNC2(gsl_ran_exponential, NoDeriv());
  TEST_FUNC(gsl_ran_exponential_pdf);
  TEST_FUNC2(gsl_cdf_exponential_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_exponential_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_exponential_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_exponential_Qinv, NoDeriv());
}

TEST_F(GSLTest, Laplace) {
  TEST_FUNC2(gsl_ran_laplace, NoDeriv());
  TEST_FUNC(gsl_ran_laplace_pdf);
  TEST_FUNC2(gsl_cdf_laplace_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_laplace_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_laplace_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_laplace_Qinv, NoDeriv());
}

TEST_F(GSLTest, ExpPow) {
  TEST_FUNC2(gsl_ran_exppow, NoDeriv());
  TEST_FUNC2(gsl_ran_exppow_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_exppow_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_exppow_Q, NoDeriv());
}

TEST_F(GSLTest, Cauchy) {
  TEST_FUNC2(gsl_ran_cauchy, NoDeriv());
  TEST_FUNC2(gsl_ran_cauchy_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_cauchy_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_cauchy_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_cauchy_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_cauchy_Qinv, NoDeriv());
}

TEST_F(GSLTest, Rayleigh) {
  TEST_FUNC2(gsl_ran_rayleigh, NoDeriv());
  TEST_FUNC2(gsl_ran_rayleigh_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_rayleigh_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_rayleigh_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_rayleigh_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_rayleigh_Qinv, NoDeriv());
}

TEST_F(GSLTest, RayleighTail) {
  TEST_FUNC2(gsl_ran_rayleigh_tail, NoDeriv());
  TEST_FUNC2(gsl_ran_rayleigh_tail_pdf, NoDeriv());
}

TEST_F(GSLTest, Landau) {
  TEST_FUNC2(gsl_ran_landau, NoDeriv());
  TEST_FUNC2(gsl_ran_landau_pdf, NoDeriv());
}

TEST_F(GSLTest, Levy) {
  TEST_FUNC2(gsl_ran_levy, NoDeriv());
  TEST_FUNC2(gsl_ran_levy_skew, NoDeriv());
}

struct GammaKnuthInfo : NoDeriv {
  bool SkipPoint(const Tuple &args) const {
    // gsl_ran_gamma_knuth doesn't work with nonpositive a.
    return args[0].number() <= 0;
  }
};

TEST_F(GSLTest, RanGamma) {
  TEST_FUNC2(gsl_ran_gamma, NoDeriv());
  TEST_FUNC2(gsl_ran_gamma_knuth, GammaKnuthInfo());
  TEST_FUNC2(gsl_ran_gamma_pdf, NoDeriv());

  Function f = GetFunction("gsl_cdf_gamma_P");
  Tuple x = MakeArgs(5, 1, 2);
  EXPECT_NEAR( 0.91791, f(x).value(), 1e-5);
  EXPECT_NEAR( 0.04104, f(x, DERIVS, BitSet("100")).deriv(), 1e-5);
  EXPECT_NEAR(-0.02052, f(x, HES, BitSet("100")).hes(), 1e-5);

  TEST_FUNC2(gsl_cdf_gamma_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_gamma_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_gamma_Qinv, NoDeriv());
}

TEST_F(GSLTest, Flat) {
  TEST_FUNC2(gsl_ran_flat, NoDeriv());
  TEST_FUNC2(gsl_ran_flat_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_flat_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_flat_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_flat_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_flat_Qinv, NoDeriv());
}

TEST_F(GSLTest, Lognormal) {
  TEST_FUNC2(gsl_ran_lognormal, NoDeriv());
  TEST_FUNC2(gsl_ran_lognormal_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_lognormal_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_lognormal_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_lognormal_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_lognormal_Qinv, NoDeriv());
}

TEST_F(GSLTest, ChiSq) {
  TEST_FUNC2(gsl_ran_chisq, NoDeriv());
  TEST_FUNC2(gsl_ran_chisq_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_chisq_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_chisq_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_chisq_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_chisq_Qinv, NoDeriv());
}

TEST_F(GSLTest, FDist) {
  TEST_FUNC2(gsl_ran_fdist, NoDeriv());
  TEST_FUNC2(gsl_ran_fdist_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_fdist_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_fdist_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_fdist_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_fdist_Qinv, NoDeriv());
}

TEST_F(GSLTest, TDist) {
  TEST_FUNC2(gsl_ran_tdist, NoDeriv());
  TEST_FUNC2(gsl_ran_tdist_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_tdist_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_tdist_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_tdist_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_tdist_Qinv, NoDeriv());
}

TEST_F(GSLTest, RanBeta) {
  // TEST_FUNC2(gsl_ran_beta, NoDeriv()); // FIXME: (GSL v2.5.0) infinity loop with arguments from POINTS
  TEST_FUNC2(gsl_ran_beta_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_beta_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_beta_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_beta_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_beta_Qinv, NoDeriv());
}

TEST_F(GSLTest, Logistic) {
  TEST_FUNC2(gsl_ran_logistic, NoDeriv());
  TEST_FUNC2(gsl_ran_logistic_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_logistic_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_logistic_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_logistic_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_logistic_Qinv, NoDeriv());
}

TEST_F(GSLTest, Pareto) {
  TEST_FUNC2(gsl_ran_pareto, NoDeriv());
  TEST_FUNC2(gsl_ran_pareto_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_pareto_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_pareto_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_pareto_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_pareto_Qinv, NoDeriv());
}

TEST_F(GSLTest, Weibull) {
  TEST_FUNC2(gsl_ran_weibull, NoDeriv());
  TEST_FUNC2(gsl_ran_weibull_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_weibull_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_weibull_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_weibull_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_weibull_Qinv, NoDeriv());
}

TEST_F(GSLTest, Gumbel1) {
  TEST_FUNC2(gsl_ran_gumbel1, NoDeriv());
  TEST_FUNC2(gsl_ran_gumbel1_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel1_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel1_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel1_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel1_Qinv, NoDeriv());
}

TEST_F(GSLTest, Gumbel2) {
  TEST_FUNC2(gsl_ran_gumbel2, NoDeriv());
  TEST_FUNC2(gsl_ran_gumbel2_pdf, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel2_P, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel2_Q, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel2_Pinv, NoDeriv());
  TEST_FUNC2(gsl_cdf_gumbel2_Qinv, NoDeriv());
}

TEST_F(GSLTest, Poisson) {
  TEST_FUNC2(gsl_ran_poisson, NoDeriv());
  TEST_FUNC2(gsl_ran_poisson_pdf, NoDeriv("k"));
  TEST_FUNC2(gsl_cdf_poisson_P, NoDeriv("k"));
  TEST_FUNC2(gsl_cdf_poisson_Q, NoDeriv("k"));
}

TEST_F(GSLTest, Bernoulli) {
  TEST_FUNC2(gsl_ran_bernoulli, NoDeriv());
  TEST_FUNC2(gsl_ran_bernoulli_pdf, NoDeriv("k"));
}

TEST_F(GSLTest, Binomial) {
  TEST_FUNC2(gsl_ran_binomial, NoDeriv("p n"));
  TEST_FUNC2(gsl_ran_binomial_pdf, NoDeriv("k p n"));
  TEST_FUNC2(gsl_cdf_binomial_P, NoDeriv("k p n"));
  TEST_FUNC2(gsl_cdf_binomial_Q, NoDeriv("k p n"));
}

struct NegativeBinomialInfo : NoDeriv {
  bool SkipPoint(const Tuple &args) const {
    // gsl_ran_negative_binomial hangs when p == 0.
    return args[0].number() == 0;
  }
};

TEST_F(GSLTest, NegativeBinomial) {
  TEST_FUNC2(gsl_ran_negative_binomial, NegativeBinomialInfo());
  TEST_FUNC2(gsl_ran_negative_binomial_pdf, NoDeriv("k"));
  TEST_FUNC2(gsl_cdf_negative_binomial_P, NoDeriv("k"));
  TEST_FUNC2(gsl_cdf_negative_binomial_Q, NoDeriv("k"));
}

struct PascalInfo : NoDeriv {
  PascalInfo() : NoDeriv("p n") {}
  bool SkipPoint(const Tuple &args) const {
    // gsl_ran_negative_binomial hangs when p == 0 and n > 0.
    return args[0].number() == 0 && args[1].number() > 0;
  }
};

TEST_F(GSLTest, Pascal) {
  TEST_FUNC2(gsl_ran_pascal, PascalInfo());
  TEST_FUNC2(gsl_ran_pascal_pdf, NoDeriv("k p n"));
  TEST_FUNC2(gsl_cdf_pascal_P, NoDeriv("k p n"));
  TEST_FUNC2(gsl_cdf_pascal_Q, NoDeriv("k p n"));
}

TEST_F(GSLTest, Geometric) {
  TEST_FUNC2(gsl_ran_geometric, NoDeriv());
  TEST_FUNC2(gsl_ran_geometric_pdf, NoDeriv("k"));
  TEST_FUNC2(gsl_cdf_geometric_P, NoDeriv("k"));
  TEST_FUNC2(gsl_cdf_geometric_Q, NoDeriv("k"));
}

TEST_F(GSLTest, HyperGeometric) {
  TEST_FUNC2(gsl_ran_hypergeometric, NoDeriv("n1 n2 t"));
  TEST_FUNC2(gsl_ran_hypergeometric_pdf, NoDeriv("k n1 n2 t"));
  TEST_FUNC2(gsl_cdf_hypergeometric_P, NoDeriv("k n1 n2 t"));
  TEST_FUNC2(gsl_cdf_hypergeometric_Q, NoDeriv("k n1 n2 t"));
}

TEST_F(GSLTest, Logarithmic) {
  TEST_FUNC2(gsl_ran_logarithmic, NoDeriv());
  TEST_FUNC2(gsl_ran_logarithmic_pdf, NoDeriv("k"));
}
}
