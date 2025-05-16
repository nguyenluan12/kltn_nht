#!/usr/bin/env bash
set -e

# lấy prefix của brew (macOS)
BREW_PREFIX=$(brew --prefix)
BOOST_INC="$BREW_PREFIX/include"
BOOST_LIB="$BREW_PREFIX/lib"

echo "🎯 Xây dựng PBWT tuần tự (C++)…"
g++-14 pbwt-src/pbwt.cpp \
  -O3 -std=c++17 \
  -I"$BOOST_INC" \
  -L"$BOOST_LIB" -lboost_iostreams -lz \
  -o pbwt

echo "🎯 Xây dựng P-PBWT song song (C++ + OpenMP)…"
g++-14 pbwt-src/p_pbwt.cpp \
  -O3 -std=c++17 -fopenmp \
  -I"$BOOST_INC" \
  -L"$BOOST_LIB" -lboost_iostreams -lz \
  -o p_pbwt

echo "🎯 Xây dựng HP-PBWT (C# .NET)…"
dotnet build hp-pbwt/source_IO_Included/HP-PBWT.csproj -c Release


echo "✅ Hoàn thành build cả 3 project."
