function j= modelJacobian_weberWrapper(t,y,growthOn,diffusionOn)
  
  [modelParams,growth]= singlecell.modelCoeffs_weber(y,growthOn,diffusionOn);
  j= singlecell.modelJacobian_weber(t,y,modelParams,growth);
end
