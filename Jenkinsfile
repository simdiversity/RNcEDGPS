pipeline {
  agent any
  stages {
    stage('Build') {
      steps {
        sh 'make build'
      }
    }
    stage('Check') {
      steps {
        sh  'make check'
      }
    }
    stage('Clean') {
      steps {
        sh 'make clean'
      }
    }
  }
}

