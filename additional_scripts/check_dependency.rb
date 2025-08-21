#!/usr/bin/env ruby

###########################################################
# Dependency Checker Script
# Checks all required dependencies for the application
###########################################################

DIR = File.dirname(__FILE__)

def print_header(title)
  puts "\n#{'=' * 50}"
  puts " #{title.upcase} ".center(50, '=')
  puts "#{'=' * 50}"
end

def print_status(message, status)
  status_color = status ? :green : :red
  puts "#{message.ljust(40)} #{status ? '[✓]'.colorize(:green) : '[✗]'.colorize(:red)}"
end

def check_ruby_version
  print_header("Ruby Version Check")
  required_version = "2.7+"
  current_version = RUBY_VERSION
  
  if Gem::Version.new(current_version) < Gem::Version.new('2.7')
    print_status("Ruby version #{current_version} (needs #{required_version})", false)
    puts "  Please install Ruby v#{required_version} using your package manager or RVM/rbenv."
    false
  else
    print_status("Ruby version #{current_version} (>= #{required_version})", true)
    true
  end
end

def check_ruby_gems
  print_header("Ruby Gems Check")
  required_gems = %w[colorize find getoptlong time fileutils tmpdir bio bio-nwk]
  missing_gems = []

  required_gems.each do |gem_name|
    installed = !Gem::Specification.find_all_by_name(gem_name).empty?
    print_status(gem_name, installed)

    unless installed
      puts "  To install: gem install #{gem_name}"
      missing_gems << gem_name
    end
  end

  if missing_gems.any?
    print "\nDo you want to install all missing gems? [Y/N]: "
    answer = STDIN.gets.strip.downcase

    if answer == "y"
      missing_gems.each do |gem_name|
        puts "Installing #{gem_name}..."
        system("gem install #{gem_name}")
      end
      puts "\n✅ All missing gems have been installed."
      true
    else
      puts "\n⚠️ Installation skipped."
      false
    end
  else
    puts "\n✅ All required gems are already installed."
    true
  end
end

def check_julia_packages
  print_header("Julia Packages Check")
  output = `julia #{DIR}/check_dependency_julia.jl --auto-install 2>&1`
  
  if $?.success? && output.strip.empty?
    print_status("All Julia packages installed", true)
    true
  else
    print_status("Some Julia packages missing", false)
    puts output
    false
  end
end

def check_r_packages
  print_header("R Packages Check")
  output = `#{DIR}/check_dependency_R.R --auto-install 2>&1`
  
  if $?.success?# && output.strip.empty?
    print_status("All R packages installed", true)
    true
  else
    print_status("Some R packages missing", false)
    puts output
    false
  end
end

def check_system_tools
  print_header("System Tools Check")
  required_tools = {
    "iqtree" => "IQ-TREE (phylogenetic software)",
    "codeml" => "CODEML (PAML package)",
    "mcmctree" => "MCMCtree (PAML package)",
    "Rscript" => "R Language",
    "julia" => "Julia Language",
    "nw_topology" => 'Newick Utilities (`conda install bioconda::newick_utils`)',
  }

  all_installed = true

  required_tools.each do |cmd, name|
    installed = system("which #{cmd} > /dev/null 2>&1")
    print_status(name, installed)
    
    unless installed
      puts "  Please install #{name}\n"
      puts "and ensure it's in your 'PATH'."
      puts "  PATH help: https://www.baeldung.com/linux/path-variable"
      all_installed = false
    end
  end

  all_installed
end



##################################################################
unless Gem::Specification.find_all_by_name('colorize').any?
  puts "Colorize gem not found. Installing..."
  success = system("gem install colorize")
  unless success
    puts "❌ Installation failed. Please install manually."
    exit 1
  end
end

require 'colorize'


# Main execution
all_checks_passed = true

all_checks_passed &= check_ruby_version
all_checks_passed &= check_ruby_gems
all_checks_passed &= check_julia_packages
all_checks_passed &= check_r_packages
all_checks_passed &= check_system_tools

if all_checks_passed
  print_header("Result")
  puts "All dependencies are properly installed!".colorize(:green)
  puts "You're ready to run the application.".colorize(:green)
else
  print_header("Result")
  puts "Some dependencies are missing or need to be updated.".colorize(:red)
  puts "Please install the missing components listed above before proceeding.".colorize(:red)
  exit 1
end
