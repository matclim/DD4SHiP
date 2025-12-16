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
  const double thinbar_x_spacing   =  x_thinbar.attr<double>("x_extra_spacing");
  const double x_offset   =  x_thinbar.attr<double>("x_offset");
  const double y_offset   =  x_thinbar.attr<double>("y_offset");
  const double extrazgap   =  x_widebar.attr<double>("extrazgap");
  const std::string calo_layer_codes = x_det.attr<std::string>("layer_codes");
  const int num_z   =  static_cast<unsigned>(calo_layer_codes.size()); 
  const int thinbar_num_x   =  x_thinbar.attr<unsigned>("num_x");

  //HPL fibre feature extraction 
  xml_dim_t    x_hplbox   = x_det.child(_Unicode(hplbox));
  
  
  //Bar definition
  Box   thinbar((x_thinbar.x()-tol)/2., (x_thinbar.y()-tol)/2.,(x_thinbar.z()-tol)/2.);
  Box   passive_layer_box((x_passive_layer.x()-tol)/2., (x_passive_layer.y()-tol)/2.,(x_passive_layer.z()-tol)/2.);
  Box   split_box((x_split.x()-tol)/2., (x_split.y()-tol)/2.,(x_split.z()-tol)/2.);
  Volume thinbar_vol("thinbar", thinbar, description.material(x_thinbar.materialStr()));
  Volume passive_layer_vol("passive_layer", passive_layer_box, description.material(x_passive_layer.materialStr()));
  Volume split_vol("split", split_box, description.material(x_split.materialStr()));
  thinbar_vol.setAttributes(description, x_thinbar.regionStr(), x_thinbar.limitsStr(), x_thinbar.visStr());
  passive_layer_vol.setAttributes(description, x_passive_layer.regionStr(), x_passive_layer.limitsStr(), x_passive_layer.visStr());
  split_vol.setAttributes(description, x_split.regionStr(), x_split.limitsStr(), x_split.visStr());


  sens.setType("calorimeter");

  // Envelope: make envelope box 'tol' bigger on each side
  Box    detbox((x_detbox.x()+tol)/2., (x_detbox.y()+tol)/2., (x_detbox.z()+tol)/2.);
  Volume detbox_vol(nam, detbox, description.air());
  detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  
//  box_vol.setVisAttributes(description.visAttributes(""));
 
  double thinlayerwidth = x_thinbar.x() * static_cast<double>(thinbar_num_x);

  Box    det_thin_layerbox((thinlayerwidth+tol)/2., (x_thinbar.y()+tol)/2., (x_thinbar.z()+tol)/2.);
  Volume det_thin_layerbox_vol("det_thin_layerbox", det_thin_layerbox, description.air());
  det_thin_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_thin_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
//  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  
  //if( x_thinbar.hasChild(_U(sensitive)) )  {
  //  sens.setType("calorimeter");
    thinbar_vol.setSensitiveDetector(sens);
  //}
  
  //if(x_hplcore.hasChild(_U(sensitive)) )  {
  //  sens.setType("calorimeter");
  //}
  //Loop for x-wise placement -> build the sensitive bar layer 
  //

  //Volume encoding
  //long DetectorCode = 9 * 1e15; 
  //long ECALCode = 1 * 1e12;

//  int DetectorCode = 9 * 1e8; 
//  int ECALCode = 1 * 1e7;

  //Thin bar layers
  int volumecode = 0;
  double xpos = -(thinlayerwidth+tol)/2.;
  for( int ix=0; ix < thinbar_num_x; ++ix )  {
    xpos += x_thinbar.x()/2.; 
    PlacedVolume pv = det_thin_layerbox_vol.placeVolume(thinbar_vol, Transform3D(rot,Position(xpos, 0e0, 0e0)));
    pv.addPhysVolID("splitcal_bar", volumecode);
    xpos += x_thinbar.x()/2. + thinbar_x_spacing; 
    volumecode++;
  }



//HPL Layers

  //Definition of layer volumes


 
  
  double z_layer = -x_detbox.z()/2.;
  Rotation3D rot_layers;


  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers

    
    std::cout << static_cast<int>(calo_layer_codes[iz]) - '0' << std::endl;
    switch(static_cast<int>(calo_layer_codes[iz]) - '0'){
    
	    case 1:{
		//Leave space for wide bars
    		z_layer += x_widebar.z();
   		break;
	    }
	    case 2:{
		//Leave space for wide bars
    		z_layer += x_widebar.z();
   		break;
	    }
	    case 3:{
		//Place thin layer vertically
    		z_layer += x_thinbar.z()/2.;
	    	rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
    	    	PlacedVolume pv_det = detbox_vol.placeVolume(det_thin_layerbox_vol, Transform3D(rot_layers,Position(x_offset,y_offset , z_layer)));
    	    	pv_det.addPhysVolID("splitcal_layer", iz);
    		z_layer += x_thinbar.z()/2.;
    		//z_layer += x_thinbar.z()+x_passive_layer.z();
   		break;
            }
	    case 4:{
		//Place thin layer horizontally	
    		z_layer += x_thinbar.z()/2.;
		rot_layers = RotationZYX(0e0, 0e0, 0e0);
    		PlacedVolume pv_det = detbox_vol.placeVolume(det_thin_layerbox_vol, Transform3D(rot_layers,Position(y_offset,x_offset, z_layer)));
        	pv_det.addPhysVolID("splitcal_layer", iz);
    		z_layer += x_thinbar.z()/2.;
    		//z_layer += x_thinbar.z()+x_passive_layer.z();
		break;
            }
	    case 5:{
		//Leave space for HPL
		z_layer += x_hplbox.z();
		break;
		   }
	    case 6:{
		//Leave space for HPL
		z_layer += x_hplbox.z();
		break;
		   }
	    case 7:{
		//Leave space for passive layers
    		z_layer += x_passive_layer.z();
    		z_layer += extrazgap;
		break;
		   }
	    case 8:{
		//Leave space for split 
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
  printout(INFO, "SplitCal Thinbars", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SplitCalThinBars,create_detector)
