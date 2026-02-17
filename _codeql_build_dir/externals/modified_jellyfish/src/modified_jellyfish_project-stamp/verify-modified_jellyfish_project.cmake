# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

if("/home/runner/work/RUFUS/RUFUS/src/modifiedJellyfish.tar.gz" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/home/runner/work/RUFUS/RUFUS/src/modifiedJellyfish.tar.gz")
  message(FATAL_ERROR "File not found: /home/runner/work/RUFUS/RUFUS/src/modifiedJellyfish.tar.gz")
endif()

if("" STREQUAL "")
  message(WARNING "File cannot be verified since no URL_HASH specified")
  return()
endif()

if("" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(VERBOSE "verifying file...
     file='/home/runner/work/RUFUS/RUFUS/src/modifiedJellyfish.tar.gz'")

file("" "/home/runner/work/RUFUS/RUFUS/src/modifiedJellyfish.tar.gz" actual_value)

if(NOT "${actual_value}" STREQUAL "")
  message(FATAL_ERROR "error:  hash of
  /home/runner/work/RUFUS/RUFUS/src/modifiedJellyfish.tar.gz
does not match expected value
  expected: ''
    actual: '${actual_value}'
")
endif()

message(VERBOSE "verifying file... done")
