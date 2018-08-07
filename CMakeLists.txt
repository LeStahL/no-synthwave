#
# Eternal Darkness - 64k intro by QM^NR4/Team210
#
# Copyright (C) 2017  Alexander Kraus <nr4@z10.info>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

cmake_minimum_required(VERSION 2.6)

project(no-synthwave)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_program(SHADER_MINIFIER NAMES shader_minifier shader_minifier.exe)
find_program(UPX NAMES upx)
add_custom_target(minification
    COMMAND ${SHADER_MINIFIER} "${CMAKE_SOURCE_DIR}/sound.frag" -o "sound.h"
    COMMAND ${SHADER_MINIFIER} "${CMAKE_SOURCE_DIR}/sound2.frag" -o "sound2.h"
    COMMAND ${SHADER_MINIFIER} "${CMAKE_SOURCE_DIR}/load210.frag" -o "load210.h"
    COMMAND ${SHADER_MINIFIER} "${CMAKE_SOURCE_DIR}/skyscrapers.frag" -o "skyscrapers.h"
    COMMAND ${SHADER_MINIFIER} "${CMAKE_SOURCE_DIR}/chipscene.frag" -o "chipscene.h"
)
add_executable(no-synthwave nswv.c)
target_link_libraries(no-synthwave X11 GL asound m pthread dl)
add_dependencies(no-synthwave minification)

add_custom_command(
    TARGET no-synthwave
    COMMAND ${UPX} --best ${PROJECT_BINARY_DIR}/no-synthwave
    POST_BUILD
)
