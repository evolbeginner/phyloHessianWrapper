#! /bin/env ruby


###########################################################
DIR = File.dirname(__FILE__)


###########################################################
$is_to_install = false


###########################################################
def err_msg(output)
  if output =~ /\w/
    puts output
    $is_to_install = true
  end
  puts
end


###########################################################
# check ruby version
puts "Checking Ruby version ......"
if RUBY_VERSION =~ /^1\. | ^2\.[0-6]/x
  puts "Your Ruby version is #{RUBY_VERSION}. Pls install Ruby v2.7+."
  $is_to_install = true
end
puts


###########################################################
# check ruby gems
puts "Checking Ruby gems ......"
required_gems = %w[parallel colorize find getoptlong time fileutils tmpdir bio bio-nwk]
missing_gems = Array.new

required_gems.each do |gem_name|
  if Gem::Specification.find_all_by_name(gem_name).empty?
    missing_gems << gem_name
    $is_to_install = true
  else
    ;
  end
end

unless missing_gems.empty?
  puts "#{missing_gems.join(',')} are NOT installed. In bash, try \'gem install #{missing_gems.join(',')}\'"
end
puts


###########################################################
# check Julia
puts "Checking Julia packages ......"
output = ` #{DIR}/check_dependency_julia.jl `
err_msg(output)


###########################################################
# check R
puts "Checking R packages ......"
output = ` #{DIR}/check_dependency_R.R `
err_msg(output)


###########################################################
# check other tools
puts "Checking other tools ......"


required_tools = {
  "iqtree"=>"iqtree",
  "nw_topology"=>"Newick Utilities",
  "codeml"=>"CODEML (PAML)",
  "mcmctree"=>"MCMCtree (PAML)",
  "Rscript"=>"R",
  "julia"=>"julia"
}

required_tools.each_pair do |k, v|
  if not system("which #{k} > /dev/null 2>&1")
    puts "#{v} NOT installed! Pls install #{v} and add its path to the environmental variable PATH (e.g., https://askubuntu.com/questions/141718/what-is-the-path-environment-variable-and-how-do-i-add-to-it)!"
    $is_to_install = true
  end
end
puts


###########################################################
unless $is_to_install
  begin
    require "colorize"
  rescue LoadError => e
    puts e.message
  ensure
    puts "Awesome! You seem to have all dependencies installed.".colorize(:red)
  end
end


