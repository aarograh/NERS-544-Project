# AUTHORS: Aaron Graham, Mike Jarrett
# PURPOSE: NERS 544 Course Project
# DATE   : April 3, 2015

# Include UseLATEX file
include(UseLATEX)

# Find Latex package
find_package(LATEX)

# Set paths
set(LATEX_SOURCE_DIR ${PROJECT_SOURCE_DIR}/doc)
set(LATEX_BINARY_DIR ${PROJECT_BINARY_DIR}/doc)

# Set source
set(LATEX_SOURCE
  introduction.tex
  conclusion.tex
)

ADD_LATEX_DOCUMENT(
  report.tex
  INPUTS ${LATEX_SOURCE}
  IMAGE_DIRS images
  DEFAULT_PDF
)

add_custom_command(TARGET pdf POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                     ${CMAKE_CURRENT_BINARY_DIR}/report.pdf
                     ${CMAKE_BINARY_DIR}/report.pdf &&
                     ${CMAKE_COMMAND} -E copy
                     ${CMAKE_CURRENT_BINARY_DIR}/report.pdf
                     ${LATEX_SOURCE_DIR}/report.pdf
)
