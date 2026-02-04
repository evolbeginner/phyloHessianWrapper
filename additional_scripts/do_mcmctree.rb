#! /usr/bin/env ruby


###############################################################
# Author: Sishuo Wang from Haiwei Luo Lab at Chinese University of Hong Kong
# E-mail: sishuowang@hotmail.ca sishuowang@cuhk.edu.hk
# Last updated: 2020-02-08
# Copyright: CC 4.0 (https://creativecommons.org/licenses/by/4.0/)
# To see the usage, run the script with '-h'

# v1.5 2024-11-08
#   rate.tre

# v1.4 2021-02-24
#   1. allows stopping before the final mcmctree (--stop_before_mcmctree)

# v1.3 2020-02-08
#   1. some bugs fixed

# v1.2 2019-12-20
#   1. allows running mcmctree for multiple trees in parallel (--tree_indir)

# v1.1 2019-11-30 (big thanks to Tianhua Liao from CUHK):
#   1. allowing running MCMCTree in parallel with "--cpu n"
#   2. disabling codeml/baseml in the first run of mcmctree


#################################################################
require 'getoptlong'
require 'parallel'
require 'time'
require 'colorize'
require 'bio-nwk'

require 'Dir'


#################################################################
$PWD = Dir.getwd

DIR ||= File.dirname(__FILE__)
FIGTREE2NWK = File.join(DIR, "figtree2tree.sh")

PAML_DIR = File.expand_path("~/software/phylo/paml/latest")
MCMCTREE_CTL = File.join(PAML_DIR, 'src', 'mcmctree.ctl')
CODEML_CTL = File.join(PAML_DIR, 'codeml.ctl')
MCMCTREE = "mcmctree"


#################################################################
seqfile = nil
tree_indir = nil
treefiles = Array.new
seqtype = nil
outdir = nil
sub_outdirs = Array.new
prefix = 'species.trees'
clock = 2
bd_paras = '1 1 0'
rgene_gamma = '1 50 1'
sigma2_gamma = '1 10 1'
alpha = 0.5
ncatG = 4
print = 2
burnin = '2000'
sampfreq = '2'
nsample = '20000'
sub_model = 'lg'
is_one_codeml = false
inBV_dir = nil
inBV_file = nil
is_inBV = false
other_args = Hash.new
cpu = 1
is_force = false
is_tolerate = false
is_stop_before_mcmctree = false

times = Array.new
$is_aa = false


#################################################################
def pf(arr)
  printf("%-50s%-80s\n", arr[0], arr[1])
end


def usage()
  STDERR.puts
  pf ['-i', "infile (in phylip format)"]
  pf ['-t', 'treefile']
  pf ['--outdir', 'outdir']
  pf ['--nucl', 'to do analysis for DNA alignments']
  pf ['--pep|prot', 'to do analysis for protein alignents']
  pf ['--clock', "the clock"]
  pf ['--BDparas|BD', "the value for BDparas (default: 1 1 0.1)"]
  pf ['--rgene|rgene_gamma', "the value for rgene_gamma (default: 1 50 1)"]
  pf ['--sigma|sigma2|sigma_gamma|sigma2_gamma', "the value for sigma2_gamma"]
  pf ['--bsn', "the value for bsn (burnin, sampfreq, nsample; default: 2000,2,20000)"]
  pf ['--alpha', "alpha in Gamma Dist"]
  pf ['--ncatG', "the no. of Gamma categories; default: 5"]
  pf ['--print', "print = ?; default: 1"]
  pf ['--nsample', "the no. of samples"]
  pf ['--sub_model', "substitution model (default: LG for protein)"]
  pf ['--one_codeml', "only a single codeml is performed."]
  pf ['--inBV', '--inbv', "specify the in.BV folder/file (automatically determined)."]
  pf ['--cpu', "cpu no."]
  pf ['--force', "remove the outdir if it exists"]
  pf ['--tolerate', "keep the outdir if it exists"]
  STDERR.puts "Please contact Sishuo Wang (sishuowang@hotmail.ca) if there are any questions. Thanks."
  exit 1
end



#################################################################
def get_ndata(seqfile)
  ndata = 0
  in_fh = File.open(seqfile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    if line =~ /^(\s+)?\d+\s+\d+$/
      ndata += 1
    end
  end
  in_fh.close
  return(ndata)
end


def sprintf_argu_line(a, b, n=14)
  str = [sprintf('%'+n.to_s+'s', a), b].join(' = ')
  return(str)
end


def prepare_paml_ctl(ctl_file, outdir, h, other_args={})
  dones = Hash.new
  lines = Array.new

  in_fh = File.open(ctl_file, 'r')
  in_fh.each_line do |line|
    line.chomp!
    if line =~ /^\s* (\w+) \s* = /x
      if h.include?($1)
        line = sprintf_argu_line($1, h[$1], 14)
        dones[$1] = ''
      end
    end
    lines << line
  end
  in_fh.close

  (h.keys - dones.keys).each do |i|
    lines << sprintf_argu_line(i, h[i], 14)
  end

  other_args.each_key do |i|
    lines << sprintf_argu_line(i, other_args[i], 14)
  end

  if Dir.exist?(outdir)
    FileUtils.cp(ctl_file, outdir)
    ctl_file = File.join(outdir, File.basename(ctl_file)) # note the ctl is changed!
  end

  out_fh = File.open(ctl_file, 'w')
  out_fh.puts lines.join("\n")
  out_fh.close
end


def run_codeml(outdir, seqtype, ndata, alpha, ncatG, cpu)
  Dir.chdir(outdir)
  `rm in.BV` if File.exist?('in.BV')

  Parallel.map(1..ndata, in_processes: cpu) do |i|
    Dir.chdir(outdir)
    c = 'tmp'+(sprintf "%04d",i)
    b = 'tmp'+(sprintf "%04d",i)+'.ctl'
    sub_outdir = File.join(outdir, c)
    `mkdir #{c}; mv #{c}.* #{c}`

    if $is_aa
      prepare_paml_ctl(File.join(sub_outdir, b), '', {'model'=>2,'aaRatefile'=>AA_RATE_FILE,'fix_alpha'=>0,'alpha'=>alpha, 'ncatG'=>ncatG, 'method'=>1})
      if ncatG == 0
        ctl_file = File.join(sub_outdir, b)
        `sed -i '/alpha/d' #{ctl_file}; sed -i '/ncatG/d' #{ctl_file}`
      end
    else
      prepare_paml_ctl(File.join(sub_outdir, b), '', {'model'=>7,'fix_alpha'=>0,'alpha'=>alpha, 'ncatG'=>ncatG}, {'method'=>1})
    end

    Dir.chdir(sub_outdir)
    if $is_aa
      `codeml #{b}`
    else
      `baseml #{b}`
    end
  end

  Dir.chdir(outdir)
  1.upto(ndata).each do |i|
    c = 'tmp'+(sprintf "%04d",i)
    `cat #{outdir}/#{c}/rst2 >> in.BV`
  end
  Dir.chdir($PWD)
end


def output_rate(print)
  output_tree = false
  k = 0
  IO.popen("grep -A1 rategram out.txt 2>/dev/null") do |grep_output|
    grep_output.each_line do |line|
      line.chomp!
      if output_tree
        File.write("rate#{k}.tre", line) # Create or overwrite the file
      end
      if line.include?("locus")
        k += 1
        output_tree = true
      else
        output_tree = false
      end
    end
  end
  
  average_branch_lengths('rate*.tre', 'rate.tre') if print == 2
  #File.symlink("rate1.tre", "rate.tre") if File.exist?("rate1.tre")
end


def average_branch_lengths(file_pattern, outfile)
  files = Dir.glob(file_pattern)
  raise "No files found matching pattern: #{file_pattern}" if files.empty?

  # Read and parse trees with error handling
  trees = files.map do |file|
    begin
      Bio::Newick.new(File.read(file)).tree
    rescue => e
      raise "Error parsing #{file}: #{e.message}"
    end
  end

  # Validate tree consistency
  n_branches = trees[0].nodes.size - 1
  trees.each_with_index do |tree, i|
    if tree.nodes.size - 1 != n_branches
      raise "Tree #{files[i]} has different topology (branch count mismatch)"
    end
  end

  # Calculate average branch lengths
  sum_branch_lengths = Array.new(n_branches+1, 0.0)
  trees.each do |tree|
    index = 0
    tree.each_edge do |node1, node2, edge|
      index += 1
      dist = edge.distance
      sum_branch_lengths[index] += dist
    end
  end

  avg_branch_lengths = sum_branch_lengths.map { |sum| sum / trees.size }

  # Create averaged tree
  avg_tree = trees.first.dup
  index = 0
  avg_tree.each_edge do |node1, node2, edge|
    index += 1
    edge.distance = avg_branch_lengths[index]
  end

  # Write output with proper file handling
  File.open(outfile, 'w') do |f|
    newick_string = avg_tree.cleanNewick
    f.write(newick_string)
  end
end



#################################################################
if __FILE__ == $0
  usage() if ARGV.empty?


  opts = GetoptLong.new(
    ['-i', GetoptLong::REQUIRED_ARGUMENT],
    ['-t', '--tree', GetoptLong::REQUIRED_ARGUMENT],
    ['--tree_indir', GetoptLong::REQUIRED_ARGUMENT],
    ['--prefix', GetoptLong::REQUIRED_ARGUMENT],
    ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--nucl', GetoptLong::NO_ARGUMENT],
    ['--pep', '--prot', GetoptLong::NO_ARGUMENT],
    ['--clock', GetoptLong::REQUIRED_ARGUMENT],
    ['--BDparas', '--BD', '--bd', GetoptLong::REQUIRED_ARGUMENT],
    ['--rgene', '--rgene_gamma', GetoptLong::REQUIRED_ARGUMENT],
    ['--sigma', '--sigma2', '--sigma_gamma', '--sigma2_gamma', GetoptLong::REQUIRED_ARGUMENT],
    ['--bsn', GetoptLong::REQUIRED_ARGUMENT],
    ['--alpha', GetoptLong::REQUIRED_ARGUMENT],
    ['--ncatG', GetoptLong::REQUIRED_ARGUMENT],
    ['--print', GetoptLong::REQUIRED_ARGUMENT],
    ['--nsample', GetoptLong::REQUIRED_ARGUMENT],
    ['--sub_model', GetoptLong::REQUIRED_ARGUMENT],
    ['--bf', GetoptLong::REQUIRED_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
    ['--force', GetoptLong::NO_ARGUMENT],
    ['--tolerate', GetoptLong::NO_ARGUMENT],
    ['--regular', GetoptLong::NO_ARGUMENT],
    ['--fastest', GetoptLong::NO_ARGUMENT],
    ['--one_codeml', GetoptLong::NO_ARGUMENT],
    ['--inBV', '--inbv', GetoptLong::REQUIRED_ARGUMENT],
    ['--stop_before_mcmctree', GetoptLong::NO_ARGUMENT],
    ['-h', GetoptLong::NO_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when '-i'
        seqfile = File.expand_path(value)
      when '-t', '--tree'
        treefiles << File.expand_path(value)
      when '--tree_indir'
        tree_indir = value
      when '--prefix'
        prefix = value
      when '--outdir'
        outdir = File.expand_path(value)
      when '--nucl'
        seqtype = 0
        $is_aa = false
        STDOUT.puts "nucl".colorize(:red)
      when '--prot', '--pep'
        seqtype = 2
        $is_aa = true
        STDOUT.puts "prot".colorize(:red)
      when '--clock'
        clock = value
        clock = 1 if clock =~ /^SR|STR$/
        clock = 2 if clock == 'IR' or clock == 'ILN'
        clock = 3 if clock == 'AR' or clock == 'GBM'
        STDOUT.puts "clock = #{clock}".colorize(:red)
      when '--BDparas', '--BD', '--bd'
        bd_paras = value.gsub(",", ' ')
      when '--rgene', '--rgene_gamma'
        rgene_gamma = value.gsub(",", ' ')
      when '--sigma', '--sigma2', '--sigma_gamma', '--sigma2_gamma'
        sigma2_gamma = value
      when '--bsn'
        arr = value.split(',')
        raise "bsn has to be x,y,z" if arr.size != 3
        burnin, sampfreq, nsample = arr
      when '--alpha'
        alpha = value.to_f
      when '--ncatG'
        ncatG = value.to_i
      when '--print'
        print = value.to_i
      when '--nsample'
        nsample = value
      when '--sub_model'
        sub_model = value.downcase
      when '--bf'
        other_args['BayesFactorBeta'] = value.to_f
      when '--cpu'
        cpu = value.to_i
      when '--force'
        is_force = true
      when '--tolerate'
        is_tolerate = true
      when '--regular'
        ;
      when '--fastest'
        clock = 1; alpha = 0.01; ncatG = 2; burnin, sampfreq, nsample = 1000, 2, 10000
        STDOUT.puts "fasttest".colorize(:blue)
      when '--one_codeml'
        is_one_codeml = true
        STDOUT.puts "one codeml!".colorize(:blue)
      when '--inBV', '--inbv'
        if File.directory?(value)
          inBV_dir = value
        else
          inBV_file = value
        end      
        is_inBV = true
        puts "is inBV".colorize(:green)
      when '--stop_before_mcmctree'
        is_stop_before_mcmctree = true
      when '-h'
        usage()
    end
  end

  puts


  #################################################################
  AA_RATE_FILE = File.join(PAML_DIR, 'dat', sub_model + '.dat')

  mkdir_with_force(outdir, is_force, is_tolerate)

  if not tree_indir.nil?
    Dir.foreach(tree_indir) do |b|
      treefiles << File.join(tree_indir, b) if b =~ /^#{prefix}/
    end
  end


  treefiles.each do |treefile|
    treefile_b = File.basename(treefile)
    outdir_suffix = treefile_b.split(/#{prefix}/)[1]
    outdir_b = ['combined', outdir_suffix].join()
    sub_outdir = File.join(outdir, outdir_b)
    sub_outdirs << sub_outdir
    mkdir_with_force(sub_outdir, true)
  end


  #################################################################
  if seqtype.nil?
    STDERR.puts "seqtype has to be given! Exiting ......"
    exit 1
  end

  ndata = get_ndata(seqfile)


  #################################################################
  outdir1 = sub_outdirs[0]
  treefile1 = treefiles[0]
  seqfile_b = File.basename(seqfile)
  treefile_b = File.basename(treefile1)

  times << Time.now # time
  puts 'mcmctree'

  sub_outdirs.each_with_index do |sub_outdir, index|
    `cp #{MCMCTREE_CTL} #{sub_outdir}`
    ctl_file = File.join(sub_outdir, File.basename(MCMCTREE_CTL))

    `cp #{seqfile} #{sub_outdir}/#{seqfile_b}`
    treefile = treefiles[index]
    `cp #{treefile} #{sub_outdir}/species.trees`

    #prepare_paml_ctl(ctl_file, sub_outdir, {'seqfile'=>seqfile_b, 'treefile'=>treefile_b, 'ndata'=>ndata, 'seqtype'=>seqtype, 'usedata'=>3, 'clock'=>clock, 'BDparas'=>bd_paras, 'rgene_gamma'=>rgene_gamma, 'sigma2_gamma'=>sigma2_gamma})
    prepare_paml_ctl(ctl_file, '', {'seqfile'=>seqfile_b, 'treefile'=>'species.trees', 'ndata'=>ndata, 'seqtype'=>seqtype, 'usedata'=>3, 'clock'=>clock, 'BDparas'=>bd_paras, 'rgene_gamma'=>rgene_gamma, 'sigma2_gamma'=>sigma2_gamma, 'burnin'=>burnin, 'sampfreq'=>sampfreq, 'nsample'=>nsample, 'alpha'=>alpha, 'ncatG'=>ncatG, 'print'=>print}, {})

    # for inBV
    if not inBV_file.nil?
      `cp #{inBV_file} #{sub_outdir}`
      next
    elsif not inBV_dir.nil?
      `cp -r #{inBV_indir}/* #{sub_outdir}`
      Dir.chdir(sub_outdir)
      `rm FigTree.tre figtree.nwk mcmctree.final`
      Dir.chdir($PWD)
      next
    end

    Dir.chdir(sub_outdir)
    `echo $$ > mcmctree.first`
    # Disable codeml
    # the last MCMCTREE might cause error "Error: ncatG" if ncatG==0 but it doesn't matter
    `mkdir disabled_bin; cd disabled_bin; touch codeml; touch baseml; chmod +x codeml; chmod +x baseml; export PATH=$PWD:$PATH; cd ../; which codeml >> #{sub_outdir}/mcmctree.first; #{MCMCTREE} mcmctree.ctl >> #{sub_outdir}/mcmctree.first`
    `mv out.BV in.BV`
    Dir.chdir($PWD)

    if is_one_codeml
      break
    end
  end

  times << Time.now; puts "Running time: " + (times[-1]-times[-2]).to_s + " seconds" # time


  #################################################################
  if $is_aa
    puts 'codeml'
  else
    puts "baseml"
  end

  if not is_inBV
    if is_one_codeml
      run_codeml(outdir1, seqtype, ndata, alpha, ncatG, cpu)
    else
      sub_outdirs.each do |sub_outdir|
        run_codeml(sub_outdir, seqtype, ndata, alpha, ncatG, cpu)
      end
    end
  end

  times << Time.now; puts "Running time: " + (times[-1]-times[-2]).to_s + " seconds" # time


  #################################################################
  if is_one_codeml
    sub_outdirs.each_with_index do |sub_outdir, index|
      next if index == 0
      `cp -r #{outdir1}/* #{sub_outdir}`
      `rm -rf #{sub_outdir}/tmp*`
      treefile = treefiles[index]
      `cp #{treefile} #{sub_outdir}/species.trees`
    end
  end


  #################################################################
  puts 'mcmctree' #start final mcmctree
  Parallel.map(sub_outdirs, in_processes: cpu) do |sub_outdir|
    Dir.chdir(sub_outdir)
    prepare_paml_ctl('mcmctree.ctl', '', {'seqfile'=>seqfile_b, 'treefile'=>'species.trees', 'ndata'=>ndata, 'seqtype'=>seqtype, 'usedata'=>"2 in.BV 3", 'clock'=>clock, 'BDparas'=>bd_paras, 'rgene_gamma'=>rgene_gamma, 'sigma2_gamma'=>sigma2_gamma, 'burnin'=>burnin, 'sampfreq'=>sampfreq, 'nsample'=>nsample, 'alpha'=>alpha, 'ncatG'=>ncatG, 'print'=>print}, other_args)
    #---------------------------------------------------------------#
    if is_stop_before_mcmctree
      puts "As required, stop before the final MCMCTree run!"; next
    end
    #---------------------------------------------------------------#

    # no matter if is_inBV or not
    `echo $$ > mcmctree.final; #{MCMCTREE} mcmctree.ctl >> mcmctree.final`
    `gzip -c mcmc.txt > mcmc.txt.gz`
    $? == 0 and `bash #{FIGTREE2NWK} -i FigTree.tre > figtree.nwk`

    #`grep rategram out.txt && grep -A1 rategram out.txt | tail -1 > rate.tre`
    output_rate(print)
    Dir.chdir($PWD)
  end


  #################################################################
  times << Time.now; puts "Running time: " + (times[-1]-times[-2]).to_s + " seconds" # time

  if $?.to_i == 0
    puts "DONE!"
  else
    puts "There was a problem. Please check it."
  end

  puts
  puts

end


