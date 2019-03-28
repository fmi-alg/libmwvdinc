#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Extended_cartesian.h>
#include <functional>

struct Dummy {
  Dummy() : value(counter++) {}
  Dummy(const Dummy & other) = default;
  Dummy(Dummy && other) = default;
  Dummy & operator=(const Dummy & other) = default;
  Dummy & operator=(Dummy && other) = default;
  int value;
  static int counter;
};

int Dummy::counter = 0;

std::ostream & operator<<(std::ostream & out, const Dummy & d) {
  out << d.value;
  return out;
}

template<typename TValue>
class Mark {
public:
  using value_type = TValue;
  using bool_ref = bool &;
  using bool_const_ref = bool const &;
public:
  Mark() : m_b(false) {}
  Mark(bool b) : m_b(b) {}
  Mark(const value_type & v) : m_b(false), m_v(v) {}
  Mark(value_type && v) : m_b(false), m_v(std::move(v)) {}
  Mark(bool b, const value_type & v) : m_b(b), m_v(v) {}
  Mark(const Mark & other) : m_b(other.m_b), m_v(other.m_v) {}
  Mark(Mark && other) : m_b(other.m_b), m_v(std::move(other.m_v)) {}
  Mark & operator=(const Mark & other) {
    m_b = other.m_b;
    m_v = other.m_v;
    return *this;
  }
  Mark & operator=(Mark && other) {
    m_b = other.m_b;
    m_v = std::move(other.m_v);
    return *this;
  }
  Mark & operator=(bool b) {
    m_b = b;
    return *this;
  }
  Mark & operator=(const value_type & v) {
    m_v = v;
    return *this;
  }
  Mark & operator=(value_type && v) {
    m_v = std::move(v);
    return *this;
  }
  // operator bool() const { return m_b; }
  // operator bool_ref() { return m_b; }
  // operator bool_const_ref() const { return m_b; }
  bool b() const { return m_b; }
  const value_type & value() const { return m_v; }
  value_type & value() { return m_v; }
public:
#define BINARY_OP(__OPNAME) \
  Mark & operator __OPNAME(bool other) { \
    m_b __OPNAME other; \
    return *this; \
  } \
  Mark & operator __OPNAME(const Mark & other) { \
    m_b __OPNAME other.b(); \
    return *this; \
  }
  BINARY_OP(&=)
  BINARY_OP(|=)
  BINARY_OP(^=)
  #undef BINARY_OP
  Mark operator!() const { return Mark(!b(), value()); }
  Mark operator~() const { return Mark(~b(), value()); }
private:
  bool m_b;
  value_type m_v;
};

#define MARK_BOOL_OPERATORS(__OPNAME) \
template<typename TValue> \
bool operator __OPNAME(bool a, const Mark<TValue> & b) { return a __OPNAME b.b(); } \
template<typename TValue> \
bool operator __OPNAME(const Mark<TValue> & a, bool b) { return a.b() __OPNAME b; } \
template<typename TValue> \
bool operator __OPNAME(const Mark<TValue> & a, const Mark<TValue> & b) { return a.b() __OPNAME b.b(); }

MARK_BOOL_OPERATORS(==)
MARK_BOOL_OPERATORS(!=)
#undef MARK_BOOL_OPERATORS

#define MARK_BINARY_BIT_OPERATORS(__OPNAME) \
template<typename TValue> \
Mark<TValue> operator __OPNAME(const Mark<TValue> & a, bool b) { return Mark<TValue>(a.b() __OPNAME b); } \
template<typename TValue> \
Mark<TValue> operator __OPNAME(bool a, const Mark<TValue> & b) { return Mark<TValue>(b.b() __OPNAME a); } \
template<typename TValue> \
Mark<TValue> operator __OPNAME(const Mark<TValue> & a, const Mark<TValue> & b) { return Mark<TValue>(a.b() __OPNAME b.b()); }

MARK_BINARY_BIT_OPERATORS(||)
MARK_BINARY_BIT_OPERATORS(&&)
MARK_BINARY_BIT_OPERATORS(|)
MARK_BINARY_BIT_OPERATORS(&)
MARK_BINARY_BIT_OPERATORS(^)
#undef MARK_BIT_OPERATORS

template<typename TValue>
std::ostream & operator<<(std::ostream & out, const Mark<TValue> & m) {
  out << "Mark{ .b=" << m.b() << "; .value()=" << m.value() << " }";
  return out;
}

//Tested with CGAL 4.12, gcc 7.3, Debian Buster

typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Epecs; //Fails with: Constructor not available for this Kernel
typedef CGAL::Extended_homogeneous<CGAL::Exact_integer>  EHEI; //Does not support roots
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec; //Does not support roots, fails with Constructor not available for this Kernel
typedef CGAL::Extended_cartesian<Epecs::RT> ECepecs; //Fails with: Segmentation fault
typedef CGAL::Simple_cartesian<Epecs::RT> SCepcs; // Fails with: Constructor not available for this Kernel
typedef Epec Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::Default_items<Kernel>::Items, Mark<Dummy>>  Nef_polyhedron;
typedef Nef_polyhedron::Plane_3  Plane_3;

int main() {
  Nef_polyhedron N1(Plane_3( 1, 0, 0,-1));
  Nef_polyhedron N2(Plane_3(-1, 0, 0,-1));
  Nef_polyhedron N3(Plane_3( 0, 1, 0,-1));
  Nef_polyhedron N4(Plane_3( 0,-1, 0,-1));
  Nef_polyhedron N5(Plane_3( 0, 0, 1,-1));
  Nef_polyhedron N6(Plane_3( 0, 0,-1,-1));
  Nef_polyhedron I1(!N1 + !N2);  // open slice in yz-plane
  Nef_polyhedron I2(N3 - !N4);   // closed slice in xz-plane
  Nef_polyhedron I3(N5 ^ N6);    // open slice in yz-plane
  Nef_polyhedron Cube1(I2 * !I1);
  Cube1 *= !I3;
  Nef_polyhedron Cube2 = N1 * N2 * N3 * N4 * N5 * N6;
  CGAL_assertion(Cube1 == Cube2);  // both are closed cube
  CGAL_assertion(Cube1 == Cube1.closure());
  CGAL_assertion(Cube1 == Cube1.regularization());
  CGAL_assertion((N1 - N1.boundary()) == N1.interior());
  CGAL_assertion(I1.closure() == I1.complement().interior().complement());
  CGAL_assertion(I1.regularization() == I1.interior().closure());

  auto c1vit = Cube1.vertices_begin();
  auto c1vend = Cube1.vertices_end();
  std::cout << "Cube vertices: " << std::endl;
  for(; c1vit != c1vend; ++c1vit) {
      std::cout << "Vertex.mark=" << c1vit->mark() << std::endl;
  }
  std::cout << "Cube halffacests: " << std::endl;
  auto c1fit = Cube1.halffacets_begin();
  auto c1fend = Cube1.halffacets_end();
  for(; c1fit != c1fend; ++c1fit) {
    std::cout << "Halffacet.mark=" << c1fit->mark() << std::endl;
  }
  return 0;
}
