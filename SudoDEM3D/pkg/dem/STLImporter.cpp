/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include"STLImporter.hpp"
#include<sudodem/pkg/dem/Shop.hpp>
#include<sudodem/lib/import/STLReader.hpp>

CREATE_LOGGER(STLImporter);

vector<shared_ptr<Body> > STLImporter::import(const char* filename)
{
	vector<Vector3r> tr;
	vector<shared_ptr<Body> > imported;

	// Load geometry
    vector<double> vtmp, ntmp; vector<int>  etmp, ftmp;
    STLReader reader; reader.tolerance=Mathr::ZERO_TOLERANCE;
    if(!reader.open(filename, back_inserter(vtmp), back_inserter(etmp), back_inserter(ftmp), back_inserter(ntmp)))
	{
		LOG_ERROR("Can't open file: " << filename);
		return imported; // zero size
	}
	for(int i=0,e=ftmp.size(); i<e; ++i)
		tr.push_back(Vector3r(vtmp[3*ftmp[i]],vtmp[3*ftmp[i]+1],vtmp[3*ftmp[i]+2]));

	// Create facets
	for(int i=0,e=tr.size(); i<e; i+=3)
	{
		Vector3r v[3]={tr[i],tr[i+1],tr[i+2]};
		Vector3r icc = Shop::inscribedCircleCenter(v[0],v[1],v[2]);
		shared_ptr<Facet> iFacet(new Facet);
		iFacet->color    = Vector3r(0.8,0.3,0.3);
		for (int j=0; j<3; ++j)
		{
				iFacet->vertices[j]=v[j]-icc;
		}
		iFacet->postLoad(*iFacet);
		shared_ptr<Body> b(new Body);
		b->state->pos=b->state->refPos=icc;
		b->state->ori=b->state->refOri=Quaternionr::Identity();
		b->shape	= iFacet;
		imported.push_back(b);
	}
	return imported;
}

