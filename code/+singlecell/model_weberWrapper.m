function dy= model_weberWrapper(t,y,N0,growthOn,diffusionOn)
  
  [modelParams,growth]= singlecell.modelCoeffs_weber(y,N0,growthOn,diffusionOn);
  dy= singlecell.model_weber(t,y,modelParams,growth);
end
