// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#ifndef INCLUDED_core_fragment_rna_RNA_MatchType_HH
#define INCLUDED_core_fragment_rna_RNA_MatchType_HH

// ObjexxFCL Headers

// you cannot #include yourself #include <protocols/rna/denovo/fragments/RNA_MatchType.hh>

// C++ Headers



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Goal: to make a fragment object that can choose fragments
// "on the fly" for RNA ab inito folding.
//
// After reading in a set of torsions from, e.g., the ribosome crystal structure,
//  should be able to generate fragments of size 1, 2, or 3, with
//  exact sequence matches, partial Y/R matches, or ignoring sequence.
//
namespace core {
namespace fragment {
namespace rna {

enum _RNA_MatchType_ { MATCH_ALL /* 0 */, MATCH_YR /* 1 */, MATCH_EXACT /* 2 */};

} //fragments
} //rna
} //protocols

#endif
