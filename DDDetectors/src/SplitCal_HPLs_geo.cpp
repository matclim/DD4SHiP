//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Climescu
// Date       : 09.11.2025
//==========================================================================
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Printout.h>
#include <iostream>
using namespace dd4hep;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  //Calo scintillator bars' feature extraction
  double       tol     = 0 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_detbox   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_widebar = x_det.child(_Unicode(widebar));
  xml_det_t    x_thinbar = x_det.child(_Unicode(thinbar));
  xml_det_t    x_passive_layer = x_det.child(_Unicode(passive_layer));
  xml_det_t    x_split = x_det.child(_Unicode(split));
  std::string  nam     = x_det.nameStr();
  //vertical bars by default
//  const double splitlayer   =  x_det.attr<int>("splitlayer");
  const double widebar_x_spacing   =  x_widebar.attr<double>("x_extra_spacing");
  const double thinbar_x_spacing   =  x_thinbar.attr<double>("x_extra_spacing");
  const double extrazgap   =  x_widebar.attr<double>("extrazgap");
  const std::string calo_layer_codes = x_det.attr<std::string>("layer_codes");
  const int num_z   =  static_cast<unsigned>(calo_layer_codes.size()); 
  const int widebar_num_x   =  x_widebar.attr<unsigned>("num_x");
  const int thinbar_num_x   =  x_thinbar.attr<unsigned>("num_x");

  //HPL fibre feature extraction 
  xml_dim_t    x_hplbox   = x_det.child(_Unicode(hplbox));
  xml_det_t    x_hplfibre = x_det.child(_Unicode(hplfibre));
  xml_det_t    x_hplcore   = x_det.child(_Unicode(hplcore));
  const double hpl_fibrethick   = x_hplfibre.thickness();
  const double hpldelta   = 2e0*x_hplfibre.rmax();
  const int    hplnum_x   = int(x_hplbox.x() / hpldelta);
  const int    hplnum_x_small = hplnum_x - 1;
//  const int    num_z   = int(2e0*x_box.z() / (delta+2*tol));
  const double hplnum_z   =  x_det.attr<int>("hpln_fibre_layers");
  
  
  //HPL definition
   
  Tube   hpl_fibre(0., x_hplfibre.rmax()-tol, (x_hplfibre.y()-tol)/2.);
  Volume hpl_fibre_vol("fibre", hpl_fibre, description.material(x_hplfibre.materialStr()));
  hpl_fibre_vol.setAttributes(description, x_hplfibre.regionStr(), x_hplfibre.limitsStr(), x_hplfibre.visStr());

  Tube   hpl_fibre_core(0., hpl_fibre.rMax()-hpl_fibrethick, (x_hplfibre.y()-tol)/2.);
  Volume hpl_fibre_core_vol("core", hpl_fibre_core, description.material(x_hplcore.materialStr()));
  hpl_fibre_core_vol.setAttributes(description, x_hplcore.regionStr(), x_hplcore.limitsStr(), x_hplcore.visStr());

  hpl_fibre_vol.placeVolume(hpl_fibre_core_vol);

  Box    hplbox((x_hplbox.x()-tol)/2., (x_hplbox.y()-tol)/2., (x_hplbox.z()-tol)/2.);
  Volume hplbox_vol(nam, hplbox, description.air());
  hplbox_vol.setAttributes(description, x_hplbox.regionStr(), x_hplbox.limitsStr(), x_hplbox.visStr());



  sens.setType("calorimeter");

  // Envelope: make envelope box 'tol' bigger on each side
  Box    detbox((x_detbox.x()+tol)/2., (x_detbox.y()+tol)/2., (x_detbox.z()+tol)/2.);
  Volume detbox_vol(nam, detbox, description.air());
  detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  

  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  

    hpl_fibre_core_vol.setSensitiveDetector(sens);

  //Volume encoding
  //long DetectorCode = 9 * 1e15; 
  //long ECALCode = 1 * 1e12;

//  int DetectorCode = 9 * 1e8; 
//  int ECALCode = 1 * 1e7;

//HPL Layers

  //Definition of layer volumes

  Box    hplbig_layer((x_hplbox.x()-tol)/2., (x_hplbox.y()-tol)/2., (x_hplfibre.rmax()-tol)/2.);
  Volume hplbig_layer_vol("splitcal_hplbig_layer", hplbig_layer, description.air());
  hplbig_layer_vol.setVisAttributes(description.visAttributes(x_hplfibre.visStr()));

  Box    hplsmall_layer((x_hplbox.x()-tol)/2., (x_hplbox.y()-tol)/2., (x_hplfibre.rmax()-tol)/2.);
  Volume hplsmall_layer_vol("splitcal_hplsmall_layer", hplsmall_layer, description.air());
  hplsmall_layer_vol.setVisAttributes(description.visAttributes(x_hplfibre.visStr()));


  int hplvolumecode =0; 
 //Build HPL layers
 
  Rotation3D hplrot(RotationZYX(0e0, 0e0, M_PI/2e0));
  for( int ix=0; ix < hplnum_x; ++ix )  {
    double x = -hplbox.x() + (double(ix)+0.5) * (hpldelta + 2e0*tol);
    PlacedVolume hplpv = hplbig_layer_vol.placeVolume(hpl_fibre_vol, Transform3D(hplrot,Position(x, 0e0, 0e0)));
    hplpv.addPhysVolID("splitcal_hplfibre", hplvolumecode);
    hplvolumecode++;
  }

  for( int ix=0; ix < hplnum_x_small; ++ix )  {
    double x = -hplbox.x() + (double(ix)+0.5) * (hpldelta + 2e0*tol) + x_hplfibre.rmax();
    PlacedVolume hplpv = hplsmall_layer_vol.placeVolume(hpl_fibre_vol, Transform3D(hplrot,Position(x, 0e0, 0e0)));
    hplpv.addPhysVolID("splitcal_hplfibre", hplvolumecode);
    hplvolumecode++;
  }

//Build the HPL Module

  for( int iz=0; iz < hplnum_z; ++iz )  {
    // leave 'tol' space between the layers
    if(iz%2 == 0){
        double z = -hplbox.z() + (double(iz)+0.5) * (2.0*tol + hpldelta);
        PlacedVolume hplpv = hplbox_vol.placeVolume(hplbig_layer_vol, Position(0e0, 0e0, z));
        hplpv.addPhysVolID("splitcal_hpl_layer", iz);
    }
    else{
        double z = -hplbox.z() + (double(iz)+0.5) * (2.0*tol + hpldelta);
        PlacedVolume hplpv = hplbox_vol.placeVolume(hplsmall_layer_vol, Position(0e0, 0e0, z));
        hplpv.addPhysVolID("splitcal_hpl_layer", iz);
    }
  }

  printout(INFO, "SHiP_HPL_Fibre_Trackers", "%s: Created %d layers of %d fibres each.", nam.c_str(), hplnum_z, hplnum_x);
  //Loop for z-wide placement -> build the calorimeter sandwich
  
  double z_layer = -x_detbox.z()/2.;
  Rotation3D rot_layers;

  int hplcount = 0;

  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers

    
    std::cout << static_cast<int>(calo_layer_codes[iz]) - '0' << std::endl;
    switch(static_cast<int>(calo_layer_codes[iz]) - '0'){
    
	    case 1:{
		//Leave space for wide bar
    		z_layer += x_widebar.z();
   		break;
	    }
	    case 2:{
		//Leave space for wide bar
    		z_layer += x_widebar.z();
   		break;
	    }
	    case 3:{
		//Leave space for thin bar
    		z_layer += x_thinbar.z();
   		break;
            }
	    case 4:{
		//Leave space for thin bar
    		z_layer += x_thinbar.z();
   		break;
            }
	    case 5:{
		//Place HPL vertically
		z_layer += x_hplbox.z()/2.;
	    	rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
    		PlacedVolume pv_det = detbox_vol.placeVolume(hplbox_vol, Transform3D(rot_layers,Position(0.,0. , z_layer)));
        	pv_det.addPhysVolID("splitcal_layer", hplcount);
    		z_layer += x_hplbox.z()/2.;
		hplcount ++;
		break;
		   }
	    case 6:{
		//Place HPL horizontally
		z_layer += x_hplbox.z()/2.;
	    	rot_layers = RotationZYX(0e0,0e0,0e0);
    		PlacedVolume pv_det = detbox_vol.placeVolume(hplbox_vol, Transform3D(rot_layers,Position(0.,0., z_layer)));
        	pv_det.addPhysVolID("splitcal_layer",hplcount);
    		z_layer += x_hplbox.z()/2.;
		hplcount ++;
		break;
		   }
	    case 7:{//Place passive layer
    		z_layer += x_passive_layer.z();
    		z_layer += extrazgap;
		break;
		   }
	    case 8:{//Place split
    		z_layer += x_split.z();
		break;
		   }
    }
    std::cout << "Zlayer Det " << z_layer << std::endl;

  } 
  DetElement   sdet  (nam, x_det.id());
  Volume       mother(description.pickMotherVolume(sdet));
  Rotation3D   rot3D (RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
  Transform3D  trafo (rot3D, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
  PlacedVolume pv = mother.placeVolume(detbox_vol, trafo);
  pv.addPhysVolID("system", x_det.id());
  sdet.setPlacement(pv);  // associate the placed volume to the detector element
  printout(INFO, "SplitCal", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SplitCalHPLs,create_detector)
