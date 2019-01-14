

FIND_PACKAGE( OpenCV REQUIRED )




ADD_EXECUTABLE( test  test.cpp  sliding_window_sample_code.cpp ring_buffer.c ＭySliding_win_SAD.cpp)
TARGET_LINK_LIBRARIES( test  ${OpenCV_LIBS} ${Boost_LIBRARIES} )



ADD_EXECUTABLE( slow_test  sad_slow_test.cpp SAD.cpp)
TARGET_LINK_LIBRARIES( slow_test  ${OpenCV_LIBS} ${Boost_LIBRARIES} )
