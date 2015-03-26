var gulp = require('gulp'),
    watch = require('gulp-watch'),
    mocha = require('gulp-mocha');

gulp.task('test', function() {
  return gulp.src('test/**/*.js', {read: false})
    .pipe(mocha({reporter: 'nyan'}));
});

gulp.task('watch', function() {
  gulp.watch('**/*.js', ['test']);
});

gulp.task('default', ['test']);
