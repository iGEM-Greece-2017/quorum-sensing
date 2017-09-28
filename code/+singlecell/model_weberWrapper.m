function dy= model_weberWrapper(t,y,growthOn,diffusionOn)
  
  [modelParams,growth]= singlecell.modelCoeffs_weber(y,growthOn,diffusionOn);
  dy= singlecell.model_weber(t,y,modelParams,growth);
end
