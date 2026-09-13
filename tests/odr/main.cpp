#include <chmath/chmath.hpp>
chm::vec3d other();
int main() {
    return chm::almost_equal(chm::length(other()), 1.0) ? 0 : 1;
}
