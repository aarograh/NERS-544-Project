# AUTHORS: Aaron Graham, Mike Jarrett
# PURPOSE: NERS 544 Course Project
# DATE   : April 3, 2015

# MACRO_ENSURE_OUT_OF_SOURCE_BUILD(<errorMessage>)

macro( MACRO_ENSURE_OUT_OF_SOURCE_BUILD _errorMessage )

string( COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" _insource )
if( _insource )
 message( SEND_ERROR "${_errorMessage}" )
 message( FATAL_ERROR
 "In-source builds are not allowed.
 It would just make a mess of the source directory.
 Please create a directory and run cmake from there, passing the path
 to this source directory as the last argument.
 This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
 Please delete them."
 )
endif( _insource )

endmacro( MACRO_ENSURE_OUT_OF_SOURCE_BUILD )
