#include <chmath/chmath.hpp>
int main() {
    const auto v = chm::normalize(chm::vec3d{3, 4, 0});
    const auto m = chm::translation(chm::vec3d{1, 2, 3});
    return chm::almost_equal(chm::length(v), 1.0) &&
                   chm::transform_point(m, chm::vec3d{}) == chm::vec3d{1, 2, 3}
               ? 0
               : 1;
}
