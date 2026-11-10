// -*- coding: utf-8 -*-
//
//  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
//  This file is part of chython.
//
//  chython is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program; if not, see <https://www.gnu.org/licenses/>.
//
// cdxconvert: read a CDX or CDXML file with Revvity's open-source ChemDraw library, write it back.
//
//   cdxconvert <input.cdx|input.cdxml> <output.cdx|output.cdxml>
//
// The input spelling is sniffed from the first bytes (VjCD0100 is CDX), the output spelling taken from
// the output extension.  Exit 0 on success; non-zero with the library's message on stderr otherwise.
// The reference reader of test_reference.py; the build recipe is in that module.

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>

#include "CDXMLParser.h"
#include "CDXStdObjects.h"

static bool ends_with(const std::string &s, const char *suffix)
{
    size_t n = std::strlen(suffix);
    return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr << "usage: cdxconvert <input.cdx|input.cdxml> <output.cdx|output.cdxml>\n";
        return 2;
    }
    std::ifstream in(argv[1], std::ios::binary);
    if (!in) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 2;
    }
    std::string data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    std::unique_ptr<CDXDocument> doc;
    try {
        if (data.compare(0, 8, "VjCD0100") == 0) {
            doc.reset(ReadCdxDiskFile(argv[1], true));
        } else {
            CDXMLParser parser;
            if (parser.XML_Parse(data.data(), static_cast<int>(data.size()), 1) == XML_STATUS_ERROR) {
                std::cerr << "CDXML parse error: " << XML_ErrorString(XML_GetErrorCode(parser))
                          << " at line " << XML_GetCurrentLineNumber(parser) << "\n";
                return 1;
            }
            doc = parser.ReleaseDocument();
        }
    } catch (const std::exception &e) {
        std::cerr << "read error: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "read error: unknown exception\n";
        return 1;
    }
    if (!doc) {
        std::cerr << "read error: no document\n";
        return 1;
    }

    std::string out(argv[2]);
    std::string message;
    bool ok = ends_with(out, ".cdx") ? WriteCdxDiskFile(doc.get(), argv[2], &message)
                                     : WriteCdxmlDiskFile(doc.get(), argv[2], &message);
    if (!ok) {
        std::cerr << "write error: " << message << "\n";
        return 1;
    }
    return 0;
}
