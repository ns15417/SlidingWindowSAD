set(OpenCV_DIR "D:/Program Files/opencv3.1/opencv/build")

FIND_PACKAGE( OpenCV REQUIRED )

ADD_EXECUTABLE( test  test.cpp  sliding_window_sample_code.cpp ring_buffer.c )
TARGET_LINK_LIBRARIES( test
                       ${OpenCV_LIBS} 
                       boost_filesystem-vc140-mt-gd-1_65_1 
                       boost_system-vc140-mt-gd-1_65_1					   
)



ADD_EXECUTABLE( slow_test  sad_slow_test.cpp SAD.cpp)
TARGET_LINK_LIBRARIES( slow_test  ${OpenCV_LIBS} ${Boost_LIBRARIES} )
