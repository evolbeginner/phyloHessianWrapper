require 'find'


################################################################################
class Dir
    def self.mkdirs(path)
        if(!File.directory?(path))
            if(!mkdirs(File.dirname(path)))
                return false;
            end
            mkdir(path)
        end
        return true
    end
end


################################################################################
def mkdir_with_force(outdir, is_force=false, is_tolerate=false)
  if outdir.class != String
    raise "outdir wrong? Exiting ......"
  end

  if ! Dir.exist?(outdir)
    `mkdir -p #{outdir}`
  else
    if is_tolerate
      ;
    elsif is_force
      `rm -rf #{outdir}`
      `mkdir -p #{outdir}`
    else
      raise "The outdir #{outdir} has already existed!"
    end
  end
end


def read_infiles(indir, suffix=nil, prefix=nil, is_all_subfolders=false)
  infiles = Array.new

  if suffix.is_a?(String)
    suffixes = [suffix]
  elsif suffix.is_a?(Array)
    suffixes = suffix
  elsif suffix.nil?
    suffixes = []
  end

  if prefix.is_a?(String)
    prefixes = [prefix]
  elsif prefix.is_a?(Array)
    prefixes = prefix
  elsif prefix.nil?
    prefixes = []
  end

  if ! is_all_subfolders
    Dir.foreach(indir) do |b|
      next if b =~ /^\./
      next if not suffixes.any?{|i| b =~ /\.#{i}$/ } unless suffixes.empty?
      next if not prefixes.any?{|i| b =~ /^#{i}/ } unless prefixes.empty?
      infiles << File.join(indir, b)
    end
  else
    Find.find(indir) do |path|
      b = File.basename(path)
      next if File.directory?(path)
      next if File.basename(path) =~ /^\./
      next if not suffixes.any?{|i| b =~ /\.#{i}$/ } unless suffixes.empty?
      next if not prefixes.any?{|i| b =~ /^#{i}/ } unless prefixes.empty?
      infiles << File.join(indir, b)
    end
  end

  return(infiles)
end


def getFilesBySuffices(indir, suffices)
  files = Array.new
  infiles = read_infiles(indir)
  infiles.each do |infile|
    if suffices.include?(File.extname(infile))
      files << infile
    end
  end
  return(files)
end


def get_file_path(file)
  path = File.symlink?(file) ? File.readlink(file) : file
  return(path)
end


