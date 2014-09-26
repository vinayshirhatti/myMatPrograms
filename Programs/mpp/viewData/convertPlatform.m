function fn = convertPlatform(fn,platform)

if strncmpi(platform,'win',3)
    fn(find(fn=='/'))='\';
elseif strncmpi(platform,'lin',3) || strncmpi(platform,'mac',3)
    fn(find(fn=='\'))='/';
end
end