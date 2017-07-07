function graphicsPause(src,event)
  if isempty(src.UserData), src.UserData= 0; end;
  if strcmp(event.Key,'space')
    if src.UserData~=1
      src.UserData= 1;
      while src.UserData~=0
        javaMethod('sleep','java.lang.Thread',100);
        drawnow;
      end
    else
      src.UserData= 0;
    end
  end
end
