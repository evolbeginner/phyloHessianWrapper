#! /usr/usr/bin/env ruby


##################################################################
#class Triez
#  require 'triez'
#  alias :include? :has_key?
#end


##################################################################
def getSpeciesNameRela(infiles)
  out2in = Hash.new
  in2out = Hash.new
  infiles.each do |infile|
    in_fh = File.open(infile, 'r')
    in_fh.each_line do |line|
      line.chomp!
      names = line.split("\t")
      outName, inName = names
      out2in[outName] = inName
      in2out[inName] = outName
    end
    in_fh.close
  end
  return([out2in, in2out])
end


def getSpeciesNameRelaWithoutFile(species_included)
  out2in = Hash.new
  in2out = Hash.new
  species_included.each_key do |taxon|
    out2in[taxon] = taxon
    in2out[taxon] = taxon
  end
  return([out2in, in2out])
end


def getGenomeBasicInfo(infile)
  h = Hash.new{|h,k|h[k]={}}
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    taxon = line_arr[0]
    h[taxon][:len] = line_arr[1].to_f
    h[taxon][:gc] = line_arr[2].to_f
  end
  in_fh.close
  return(h)
end


##################################################################
def read_list(infile=nil, field=1, sep="\t", is_arr=false)
  genes = is_arr ? Array.new : Hash.new
  return(genes) if infile.nil?

  in_fh = infile == '-' ? STDIN : File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split(sep)
    item = line_arr[field-1]
    # get_the_next_one
    next_item = line_arr[line_arr.size-field]
    genes.is_a?(Array) ? genes << item : genes[item] = next_item
  end
  in_fh.close if infile != '-'
  return(genes)
end


def read_lines(infile=nil)
  lines = Array.new
  return(lines) if infile.nil?
  in_fh = infile == '-' ? STDIN : File.open(infile, 'r')
  lines = in_fh.readlines.map!{|i|i.chomp}
  in_fh.close if infile != '-'
  return(lines)
end


def readTbl(infile, fields=[1, 2], sep="\t")
  field1, field2 = fields
  a2b, b2a = Hash.new, Hash.new
  in_fh = infile == '-' ? STDIN : File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split(sep)
    field1 = 1 and field2 = 1 if line_arr.size < 2 # to deal w/ fam.list w/ a single line
    items = line_arr.values_at(field1-1, field2-1)
    a2b[items[0]] = items[1]
    b2a[items[1]] = items[0]
  end
  in_fh.close if infile != '-'
  return([a2b, b2a])
end

alias read_tbl readTbl


def read_list_2_arr(infile=nil, field=1, sep="\t")
  genes = Array.new
  return(genes) if infile.nil?
  in_fh = infile == '-' ? STDIN : File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split(sep)
    item = line_arr[field-1]
    genes << item
  end
  in_fh.close if infile != '-'
  return(genes)
end


def reformat_wrapped(s, width=78, num_of_space=0)
  # https://www.safaribooksonline.com/library/view/ruby-cookbook/0596523696/ch01s15.html
  lines = []
  line = ""
  s.split(/\s+/).each do |word|
    if line.size + word.size >= width
      lines << line
      line = word
    elsif line.empty?
      line = word
    else
      line << " " << word
    end
  end
  lines << line if line

  lines.each_with_index do |line, index|
    break if num_of_space == 0
    next if index == 0
    lines[index] = ' ' * num_of_space + line
  end

  return lines.join "\n"
end


def processbar(count, total)
  maxlen = 80
  prop = (count.to_f/total)*100
  jing = '#' * (prop*maxlen/100)
  printf "\r%-80s%s", jing, prop.to_s+'%'
end


def getCorename(infile, is_strict=false)
  b = File.basename(infile)
  if not is_strict
    if b =~ /\./
      b =~ /(.+)\..+$/
    else
      b =~ /^(.+)$/
    end
  else
    b =~ /^([^.]+)/
  end
  c = $1
  return(c)
end


def cp2dir(in_items, out)
  require 'fileutils'
  [in_items].flatten.each do |i|
    FileUtils.cp(i, out)
  end
end


def errMsg(msg, color=:red)
  require 'colorize'
  STDERR.puts msg.colorize(color)
  exit 1
end


