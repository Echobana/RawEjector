from cx_Freeze import setup, Executable

executables = [Executable('EGPlot.py')]

excludes = ['unicodedata', 'logging', 'unittest', 'email', 'html', 'http', 'urllib',
            'xml', 'bz2']

zip_include_packages = ['collections', 'encodings', 'importlib', 'unicodedata']

includes = ['unicodedata']

options = {
    'build_exe': {
        'include_msvcr': True,
        'excludes': excludes,
    }
}

setup(name='EGPlot',
      version='0.0.3',
      description='My Hello World App!',
      executables=executables,
      options=options)
