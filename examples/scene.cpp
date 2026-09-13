#include <chmath/geometry/frustum.hpp>
#include <chmath/geometry/intersection.hpp>
#include <cstdio>
int main() {
    using namespace chm;
    const auto q = from_axis_angle(vec3f{0, 1, 0}, radians(30.0f));
    if (!q)
        return 1;
    const auto model = compose(trs<float>{{0, 0, -5}, *q, {1, 1, 1}});
    const auto view = look_at(vec3f{0, 2, 5}, vec3f{0, 0, -5}, vec3f{0, 1, 0});
    const auto projection = perspective(radians(60.0f), 16.0f / 9.0f, 0.1f, 1000.0f);
    if (!view || !projection)
        return 2;
    const auto vp = *projection * to_matrix(*view);
    const auto fr = extract_frustum(vp);
    if (!fr)
        return 3;
    const auto bounds = transform_bounds(model, aabb3f{{-1, -1, -1}, {1, 1, 1}});
    const auto screen =
        project(transform_point(model, vec3f{}), vp, viewport<float>{0, 0, 1920, 1080});
    if (!screen)
        return 4;
    std::printf("center=(%.2f, %.2f), visible=%s, vec3=%zu bytes, affine=%zu bytes\n",
                static_cast<double>((*screen)[0]), static_cast<double>((*screen)[1]),
                classify(*fr, bounds) != containment::outside ? "yes" : "no", sizeof(vec3f),
                sizeof(affine3f));
    return 0;
}
