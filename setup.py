from numpy.distutils.core import Extension

pyenkf = Extension(name                   = 'pyenkf',
                   sources                = ['src/pyenkf/enksrf.F90'],
                   extra_f90_compile_args = ["-O3","-fdefault-real-8"],
                   extra_link_args        = ["-lopenblas"],
                   )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name         = 'pyenkf',
          version      = "0.0.1",
          description  = "Ensemble Kalman Square Root Filter API",
          author       = "J.J.D. Hooghiem",
          author_email = "joramjd@gmail.com",
          url          = "https://github.com/JJDHooghiem/pyenkf",
          package_dir  = {"":"src"},
          ext_modules  = [pyenkf]
          )
