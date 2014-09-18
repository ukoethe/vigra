VIGRA_ADD_TEST(test_multidistance test.cxx LIBRARIES vigraimpex)

VIGRA_COPY_TEST_DATA(
    blatt.xv
    raw_skeleton.xv
    skeleton_length.xv
    skeleton_length_greater_100.xv
    skeleton_length_greater_50_percent.xv
    skeleton_salience.xv
    skeleton_salience_greater_20.xv
    skeleton_salience_greater_60_percent.xv
    skeleton_topology.xv
    skeleton_topology_without_center.xv
)
