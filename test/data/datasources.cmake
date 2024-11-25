# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# This includes `cmake/test/declare_datasource.cmake`, which provides the `declare_datasource` function.
include (test/declare_datasource)

# Makes all files in this directory and subdirectories available to the build, i.e., copying them to <build>/data/.
# `datasources.cmake`, `README.md`, and files ending in `.license` are ignored.
# You may organise your data in subdirectories, but each file must have a unique name.
include (test/add_local_data)

# Downloads the file to <build>/data/downloaded.fasta
# The checksum of the downloaded file is checked to ensure that:
# * The download was not corrupted.
# * The file was not tampered with.
#
# Some ways to obtain the checksum:
#
#   wget --quiet -O- https://ftp.seqan.de/app-template/downloaded.fasta | sha256sum
#   c3cb990ca1a25c7e31be3c6c2d009238d9ac9a44b2b7c143753c1e2881699077  -
#   ^------------------------- checksum ---------------------------^  ^ file (- is stdin, the file was piped)
#
# You can also use curl:
#   curl --silent https://ftp.seqan.de/app-template/downloaded.fasta | sha256sum
#
# If you have the file locally, you can use `sha256sum` directly:
#   sha256sum downloaded.fasta
#   c3cb990ca1a25c7e31be3c6c2d009238d9ac9a44b2b7c143753c1e2881699077  downloaded.fasta
#declare_datasource (FILE downloaded.fasta # This is a custom name. It does not have to match the file name.
#                    URL https://ftp.seqan.de/app-template/downloaded.fasta
#                    URL_HASH SHA256=c3cb990ca1a25c7e31be3c6c2d009238d9ac9a44b2b7c143753c1e2881699077
#)

#declare_datasource (FILE meta.bin # This is a custom name. It does not have to match the file name.
#                    URL ${CMAKE_SOURCE_DIR}/test/data/meta.bin
#                    URL_HASH SHA256=38a45a3aad8f824e8f19c4971aba53e138d5075e3c53686b402b525ab720c546
#)

#declare_datasource (FILE truth.gff 
#                    URL ${CMAKE_SOURCE_DIR}/test/data/truth.gff
#                    URL_HASH SHA256=0d92beee9e6646f61913b0d606f4e10b924b0697d98c71cf15af441d8992f65a
#)

#declare_datasource (FILE test.gff 
#                    URL ${CMAKE_SOURCE_DIR}/test/data/test.gff
#                    URL_HASH SHA256=6e227ff7f1ef8a387c3175d65214c5c1ce64c5cd5d4615fb48c6777f844e925e
#)